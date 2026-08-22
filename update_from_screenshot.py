"""Read the fixed ITRTG event grid from a screenshot and update an optimizer.

The seasonal artwork and item names change, but the 3-by-4 generator grid keeps
the same normalized screen coordinates.  Each cell is OCR'd independently so
text from adjacent columns cannot be interleaved by the Windows OCR engine.
"""

from __future__ import annotations

import argparse
import asyncio
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
import os
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

try:
    from PIL import Image, ImageChops, ImageEnhance, ImageFilter
except ImportError:
    Image = None

OCR_IMPORT_ERROR: ImportError | None = None
try:
    # Maintained, split PyWinRT packages.
    import winrt.windows.foundation
    import winrt.windows.globalization
    import winrt.windows.storage.streams
    from winrt.windows.graphics.imaging import BitmapDecoder
    from winrt.windows.media.ocr import OcrEngine
    from winrt.windows.storage import FileAccessMode, StorageFile
except ImportError as modern_error:
    try:
        # Compatibility with older installations of the monolithic package.
        from winsdk.windows.graphics.imaging import BitmapDecoder
        from winsdk.windows.media.ocr import OcrEngine
        from winsdk.windows.storage import FileAccessMode, StorageFile
    except ImportError:
        BitmapDecoder = OcrEngine = FileAccessMode = StorageFile = None
        OCR_IMPORT_ERROR = modern_error


ITEM_COUNT = 12

# Normalized coordinates measured from multiple 2048x1188/1189 screenshots.
# Only the right-hand text portion of each item cell is included, avoiding most
# seasonal artwork while retaining the name, level, speed, and count.
COLUMN_RANGES = ((0.384, 0.524), (0.592, 0.731), (0.800, 0.939))
ROW_RANGES = ((0.367, 0.500), (0.518, 0.651), (0.670, 0.802), (0.821, 0.954))
CELL_STRIPS = {
    # (left, right, top, bottom), relative to the text portion of a cell.
    "name": (0.00, 1.00, 0.00, 0.27),
    # Numeric strips exclude the labels. Windows OCR otherwise sometimes sees
    # "Level:" but treats the visually separated value as another region.
    "level": (0.28, 0.76, 0.22, 0.47),
    "speed": (0.28, 0.76, 0.43, 0.69),
    "count": (0.28, 1.00, 0.66, 0.93),
}

LABEL_PATTERNS = {
    "level": re.compile(r"\bLeve[lI1]\s*[:;.]?\s*", re.IGNORECASE),
    "speed": re.compile(r"\bSpee[dcl]\s*[:;.]?\s*", re.IGNORECASE),
    "count": re.compile(r"\bCou[nm][tI1]\s*[:;.]?\s*", re.IGNORECASE),
}
NUMBER_PATTERN = re.compile(
    r"([0-9][0-9,]*(?:\.[0-9]+)?)\s*"
    r"(thousand|million|billion|trillion)?",
    re.IGNORECASE,
)
MAGNITUDES = {
    None: 1,
    "thousand": 1_000,
    "million": 1_000_000,
    "billion": 1_000_000_000,
    "trillion": 1_000_000_000_000,
}

@dataclass(frozen=True)
class ItemStats:
    name: str
    level: int
    speed: int
    count: int
    ocr_text: str


@dataclass(frozen=True)
class ScreenshotData:
    items: list[ItemStats]
    time_left: tuple[int, int, int, int] | None
    full_ocr_text: str


def require_ocr_dependencies() -> None:
    if Image is None:
        raise RuntimeError(
            "Pillow is not installed. Run: python -m pip install -r requirements-ocr.txt"
        )
    if OCR_IMPORT_ERROR is not None:
        raise RuntimeError(
            "Windows OCR bindings are not installed. Run: "
            "python -m pip install -r requirements-ocr.txt"
        ) from OCR_IMPORT_ERROR


async def recognize_file(file_path: Path) -> str:
    file = await StorageFile.get_file_from_path_async(str(file_path.resolve()))
    stream = await file.open_async(FileAccessMode.READ)
    try:
        decoder = await BitmapDecoder.create_async(stream)
        software_bitmap = await decoder.get_software_bitmap_async()
        engine = OcrEngine.try_create_from_user_profile_languages()
        if not engine:
            raise RuntimeError("Windows could not create an OCR engine for your language profile.")
        result = await engine.recognize_async(software_bitmap)
        return result.text
    finally:
        close = getattr(stream, "close", None)
        if close:
            close()


def normalized_box(
    width: int, height: int, x_range: tuple[float, float], y_range: tuple[float, float]
) -> tuple[int, int, int, int]:
    return (
        round(width * x_range[0]),
        round(height * y_range[0]),
        round(width * x_range[1]),
        round(height * y_range[1]),
    )


def prepare_crop(image, box: tuple[int, int, int, int], monochrome: bool = False):
    crop = image.crop(box).convert("RGB")
    crop = crop.resize((crop.width * 3, crop.height * 3), Image.Resampling.LANCZOS)
    if monochrome:
        # The UI text is near-white. Inverting a luminance mask gives Windows
        # OCR conventional black text on white and removes most event artwork.
        gray = crop.convert("L")
        return gray.point(lambda value: 0 if value >= 175 else 255).convert("RGB")
    crop = ImageEnhance.Contrast(crop).enhance(1.15)
    return crop.filter(ImageFilter.SHARPEN)


def trim_white(crop):
    difference = ImageChops.difference(crop, Image.new("RGB", crop.size, "white"))
    bounds = difference.getbbox()
    return crop.crop(bounds) if bounds else crop


def join_label_and_value(label_crop, value_crop):
    """Remove the UI's wide label/value gap to help Windows OCR group them."""
    label_crop = trim_white(label_crop)
    value_crop = trim_white(value_crop)
    padding = 15
    height = max(label_crop.height, value_crop.height) + padding * 2
    joined = Image.new(
        "RGB", (label_crop.width + value_crop.width + padding * 3, height), "white"
    )
    joined.paste(label_crop, (padding, padding))
    joined.paste(value_crop, (label_crop.width + padding * 2, padding))
    return joined


def parse_number(value: str) -> int:
    # OCR occasionally inserts a space inside a number over detailed artwork.
    value = re.sub(r"(?<=[\d,.])\s+(?=[\d,.])", "", value)
    match = NUMBER_PATTERN.search(value)
    if not match:
        # In isolated numeric fields, Windows OCR commonly renders 1 as I/l
        # and 0 as O. Require a token boundary so label text cannot become a
        # fabricated number.
        ambiguous = re.search(r"(?<![A-Za-z])([Il!Oo](?:\s*[Il!Oo])*)(?![A-Za-z])", value)
        if not ambiguous:
            raise ValueError(f"no number found in {value!r}")
        normalized = re.sub(r"\s+", "", ambiguous.group(1)).translate(
            str.maketrans({"I": "1", "l": "1", "!": "1", "O": "0", "o": "0"})
        )
        return int(normalized)
    number = float(match.group(1).replace(",", ""))
    magnitude = match.group(2).lower() if match.group(2) else None
    return round(number * MAGNITUDES[magnitude])


def parse_labeled_number(text: str, label: str) -> int:
    match = LABEL_PATTERNS[label].search(text)
    if not match:
        raise ValueError(f"missing {label} label")
    return parse_number(text[match.end() :])


def parse_name(text: str) -> str:
    first_label = min(
        (match.start() for pattern in LABEL_PATTERNS.values() if (match := pattern.search(text))),
        default=len(text),
    )
    name = " ".join(text[:first_label].split()).strip(" -:;|.")
    if not name or not re.search(r"[A-Za-z]", name):
        raise ValueError("missing item name")
    return name


def parse_strip_number(text: str, label: str) -> int:
    """Parse a stat strip, tolerating a damaged or omitted label."""
    try:
        return parse_labeled_number(text, label)
    except ValueError:
        return parse_number(text)


def first_parsed(candidates: list[str], parser, description: str):
    errors = []
    for candidate in candidates:
        try:
            return parser(candidate)
        except ValueError as error:
            errors.append(str(error))
    raise ValueError(f"{description}: {'; '.join(errors)}")


def majority_parsed(candidates: list[str], parser):
    values = []
    for candidate in candidates:
        try:
            values.append(parser(candidate))
        except ValueError:
            pass
    if not values:
        return None
    ranking = Counter(values).most_common()
    if len(ranking) == 1 or ranking[0][1] > ranking[1][1]:
        return ranking[0][0]
    return None


def consensus_parsed(candidates: list[str], parser):
    values = set()
    for candidate in candidates:
        try:
            values.add(parser(candidate))
        except ValueError:
            pass
    if len(values) == 1:
        return values.pop()
    return None


def glyph_components(image) -> list:
    """Return left-to-right digit-sized connected components from a binary crop."""
    gray = image.convert("L")
    width, height = gray.size
    pixels = gray.load()
    dark = {(x, y) for y in range(height) for x in range(width) if pixels[x, y] < 128}
    components = []
    while dark:
        seed = dark.pop()
        stack = [seed]
        points = [seed]
        while stack:
            x, y = stack.pop()
            for neighbor in (
                (x - 1, y - 1), (x, y - 1), (x + 1, y - 1),
                (x - 1, y),                     (x + 1, y),
                (x - 1, y + 1), (x, y + 1), (x + 1, y + 1),
            ):
                if neighbor in dark:
                    dark.remove(neighbor)
                    stack.append(neighbor)
                    points.append(neighbor)
        xs = [point[0] for point in points]
        ys = [point[1] for point in points]
        box = (min(xs), min(ys), max(xs) + 1, max(ys) + 1)
        if box[3] - box[1] >= height * 0.35 and len(points) >= 40:
            components.append((box[0], box[2], image.crop(box)))
    components.sort(key=lambda value: value[0])
    # Bright seasonal artwork can survive thresholding near the far edge of a
    # strip. Numeric glyphs form one tight horizontal cluster, so discard
    # components following a large gap.
    clustered = []
    previous_right = None
    for component in components:
        left, right, glyph = component
        if previous_right is not None and left - previous_right > height * 0.55:
            break
        clustered.append(glyph)
        previous_right = right
    return clustered


def normalize_glyph(glyph):
    canvas = Image.new("L", (40, 60), 255)
    glyph = glyph.convert("L")
    glyph.thumbnail((32, 52), Image.Resampling.NEAREST)
    canvas.paste(glyph, ((40 - glyph.width) // 2, (60 - glyph.height) // 2))
    return canvas


def glyph_distance(left, right) -> float:
    left_image = normalize_glyph(left)
    right_image = normalize_glyph(right)
    left_getter = getattr(left_image, "get_flattened_data", left_image.getdata)
    right_getter = getattr(right_image, "get_flattened_data", right_image.getdata)
    left_pixels = list(left_getter())
    right_pixels = list(right_getter())
    different = sum((a < 128) != (b < 128) for a, b in zip(left_pixels, right_pixels))
    ink = sum((a < 128) or (b < 128) for a, b in zip(left_pixels, right_pixels))
    return different / max(ink, 1)


def recognize_from_templates(image, templates: dict[str, list]) -> int:
    components = glyph_components(image)
    if not components:
        raise ValueError("no digit-shaped components found")
    digits = []
    for component_index, component in enumerate(components, start=1):
        scores = sorted(
            (
                min(glyph_distance(component, example) for example in examples),
                digit,
            )
            for digit, examples in templates.items()
            if examples
        )
        if not scores or scores[0][0] > 0.60:
            best = f"{scores[0][1]} at {scores[0][0]:.3f}" if scores else "none"
            raise ValueError(
                f"digit {component_index} did not match a learned screenshot template "
                f"(best: {best})"
            )
        if len(scores) > 1 and scores[1][0] - scores[0][0] < 0.04:
            raise ValueError(
                f"digit {component_index} template match was ambiguous "
                f"({scores[0][1]}={scores[0][0]:.3f}, "
                f"{scores[1][1]}={scores[1][0]:.3f})"
            )
        digits.append(scores[0][1])
    return int("".join(digits))


def parse_time_left(text: str) -> tuple[int, int, int, int] | None:
    match = re.search(
        r"Event\s+time\s+left\s*[:;.]\s*"
        r"(?:(\d+)\s*days?\s*[,]?\s*)?"
        r"(\d+)\s*:\s*(\d+)\s*:\s*(\d+)",
        text,
        re.IGNORECASE,
    )
    if not match:
        return None
    days, hours, minutes, seconds = match.groups()
    return int(days or 0), int(hours), int(minutes), int(seconds)


async def extract_screenshot(file_path: Path) -> ScreenshotData:
    require_ocr_dependencies()
    if not file_path.is_file():
        raise FileNotFoundError(f"Could not find screenshot: {file_path}")

    with Image.open(file_path) as source:
        source.load()
        full_text = await recognize_file(file_path)
        cell_records = []
        failures: list[str] = []

        with TemporaryDirectory(prefix="itrtg-ocr-") as temp_dir:
            for row, y_range in enumerate(ROW_RANGES):
                for column, x_range in enumerate(COLUMN_RANGES):
                    index = row * 3 + column
                    box = normalized_box(source.width, source.height, x_range, y_range)
                    left, top, right, bottom = box
                    cell_height = bottom - top
                    parts: dict[str, list[str]] = {}
                    numeric_images = {}
                    cell_width = right - left
                    for part, (x_start, x_end, y_start, y_end) in CELL_STRIPS.items():
                        strip_box = (
                            round(left + cell_width * x_start),
                            round(top + cell_height * y_start),
                            round(left + cell_width * x_end),
                            round(top + cell_height * y_end),
                        )
                        parts[part] = []
                        for variant, monochrome in (("color", False), ("mono", True)):
                            crop = prepare_crop(source, strip_box, monochrome=monochrome)
                            crop_path = (
                                Path(temp_dir) / f"cell-{index + 1}-{part}-{variant}.png"
                            )
                            crop.save(crop_path)
                            parts[part].append(await recognize_file(crop_path))
                            if monochrome and part in {"level", "speed", "count"}:
                                numeric_images[part] = crop.copy()
                        if part in {"level", "speed"}:
                            label_box = (
                                left,
                                strip_box[1],
                                round(left + cell_width * 0.31),
                                strip_box[3],
                            )
                            for variant, monochrome in (("color", False), ("mono", True)):
                                joined = join_label_and_value(
                                    prepare_crop(source, label_box, monochrome=monochrome),
                                    prepare_crop(source, strip_box, monochrome=monochrome),
                                )
                                joined_path = (
                                    Path(temp_dir)
                                    / f"cell-{index + 1}-{part}-joined-{variant}.png"
                                )
                                joined.save(joined_path)
                                parts[part].append(await recognize_file(joined_path))
                    try:
                        name = first_parsed(parts["name"], parse_name, "missing item name")
                    except ValueError as error:
                        failures.append(
                            f"cell {index + 1} (row {row + 1}, column {column + 1}): "
                            f"{error}; OCR={parts!r}"
                        )
                        continue
                    cell_records.append(
                        {
                            "index": index,
                            "name": name,
                            "level": majority_parsed(
                                parts["level"],
                                lambda text: parse_strip_number(text, "level"),
                            ),
                            "speed": majority_parsed(
                                parts["speed"],
                                lambda text: parse_strip_number(text, "speed"),
                            ),
                            # Counts must agree across preprocessing variants.
                            # A disagreement is resolved below using glyphs
                            # learned from unambiguous fields in this image.
                            "count": consensus_parsed(
                                parts["count"],
                                lambda text: parse_strip_number(text, "count"),
                            ),
                            "parts": parts,
                            "numeric_images": numeric_images,
                        }
                    )

    if failures:
        raise ValueError("Could not reliably read the generator grid:\n  " + "\n  ".join(failures))
    if len(cell_records) != ITEM_COUNT:
        raise ValueError(f"Expected {ITEM_COUNT} item cells, read {len(cell_records)}")

    templates: dict[str, list] = {str(digit): [] for digit in range(10)}
    for record in cell_records:
        for field in ("level", "speed"):
            value = record[field]
            components = glyph_components(record["numeric_images"][field])
            # Reject syntactically valid partial OCR (for example, 13 read as
            # a lone "I" -> 1) when the image visibly contains more glyphs.
            if value is not None and len(components) != len(str(value)):
                record[field] = value = None
            if value is not None and len(components) == len(str(value)):
                for digit, component in zip(str(value), components):
                    templates[digit].append(component)
        count_has_magnitude = any(
            re.search(r"thousand|million|billion|trillion", candidate, re.IGNORECASE)
            for candidate in record["parts"]["count"]
        )
        if record["count"] is not None and not count_has_magnitude:
            components = glyph_components(record["numeric_images"]["count"])
            if len(components) == len(str(record["count"])):
                for digit, component in zip(str(record["count"]), components):
                    templates[digit].append(component)

    items: list[ItemStats] = []
    template_failures = []
    for record in cell_records:
        for field in ("level", "speed", "count"):
            if record[field] is None:
                try:
                    record[field] = recognize_from_templates(
                        record["numeric_images"][field], templates
                    )
                except ValueError as error:
                    template_failures.append(
                        f"cell {record['index'] + 1} {field}: {error}; "
                        f"OCR={record['parts'][field]!r}"
                    )
        if record["level"] is not None and record["speed"] is not None:
            items.append(
                ItemStats(
                    name=record["name"],
                    level=record["level"],
                    speed=record["speed"],
                    count=record["count"],
                    ocr_text=" | ".join(
                        f"{key}={value!r}" for key, value in record["parts"].items()
                    ),
                )
            )

    if template_failures or len(items) != ITEM_COUNT:
        raise ValueError(
            "Could not reliably read the generator grid:\n  " + "\n  ".join(template_failures)
        )

    return ScreenshotData(items=items, time_left=parse_time_left(full_text), full_ocr_text=full_text)


def format_time(time_left: tuple[int, int, int, int] | None) -> str:
    if time_left is None:
        return "not recognized"
    days, hours, minutes, seconds = time_left
    return f"{days} days, {hours:02}:{minutes:02}:{seconds:02}"


def print_data(data: ScreenshotData, show_ocr: bool = False) -> None:
    print(f"Event time left: {format_time(data.time_left)}")
    print()
    print(f"{'#':>2}  {'Item':<20} {'Level':>5} {'Speed':>5} {'Count':>16}")
    print("--  -------------------- ----- ----- ----------------")
    for index, item in enumerate(data.items, start=1):
        print(f"{index:>2}  {item.name:<20.20} {item.level:>5} {item.speed:>5} {item.count:>16,}")
    if show_ocr:
        print("\n--- Full-screen OCR ---")
        print(data.full_ocr_text)
        print("\n--- Cell OCR ---")
        for index, item in enumerate(data.items, start=1):
            print(f"{index:>2}: {item.ocr_text!r}")


def normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def cpp_resource_names(content: str) -> list[str]:
    match = re.search(r"resourceNames\s*=\s*\{(.*?)\};", content, re.DOTALL)
    if not match:
        return []
    return re.findall(r'"([^"]+)"', match.group(1))


def validate_target_event(content: str, items: list[ItemStats]) -> None:
    expected = cpp_resource_names(content)
    if len(expected) < 6:
        raise ValueError("Could not read the target C++ resourceNames array.")
    mismatches = [
        (index + 1, expected[index], items[index].name)
        for index in range(6)
        if normalize_name(expected[index]) != normalize_name(items[index].name)
    ]
    if mismatches:
        details = "; ".join(
            f"cell {index}: target={wanted!r}, screenshot={found!r}"
            for index, wanted, found in mismatches
        )
        raise ValueError(f"Screenshot does not match the target event ({details}).")


def replace_once(content: str, pattern: str, replacement: str, description: str) -> str:
    updated, count = re.subn(pattern, replacement, content, count=1, flags=re.DOTALL)
    if count != 1:
        raise ValueError(f"Could not uniquely locate {description} in the target C++ file.")
    return updated


def update_cpp(cpp_path: Path, data: ScreenshotData, force: bool = False) -> None:
    if not cpp_path.is_file():
        raise FileNotFoundError(f"Could not find C++ target: {cpp_path}")
    content = cpp_path.read_text(encoding="utf-8")
    if not force:
        validate_target_event(content, data.items)

    levels = [item.level for item in data.items]
    speeds = [item.speed for item in data.items]
    counts = [item.count for item in data.items]

    levels_text = f"""array<int, 25> currentLevels = {{
    // Current Production Levels
    {levels[0]}, {levels[1]}, {levels[2]},
    {levels[3]}, {levels[4]}, {levels[5]},
    {levels[6]}, {levels[7]}, {levels[8]},
    {levels[9]}, {levels[10]}, {levels[11]}, 
    
    // Current Speed Levels
    {speeds[0]}, {speeds[1]}, {speeds[2]},
    {speeds[3]}, {speeds[4]}, {speeds[5]},
    {speeds[6]}, {speeds[7]}, {speeds[8]},
    {speeds[9]}, {speeds[10]}, {speeds[11]},
    
    0 // Dummy placeholder. Keep 0
}};"""

    counts_text = f"""array<double, 12> resourceCounts = {{
    // Current resource counts
    {counts[0]}, {counts[1]}, {counts[2]}, 
    {counts[3]}, {counts[4]}, {counts[5]}, 
    {counts[6]}/((500.0+DLs)/5.0), {counts[7]}, {counts[8]},
    {counts[9]}/(0.5+AL/100.0), {counts[10]}, {counts[11]}/(UNLOCKED_PETS/100.0)
}};"""

    updated = replace_once(
        content,
        r"array<int,\s*25>\s+currentLevels\s*=\s*\{.*?\};",
        levels_text,
        "currentLevels",
    )
    updated = replace_once(
        updated,
        r"array<double,\s*12>\s+resourceCounts\s*=\s*\{.*?\};",
        counts_text,
        "resourceCounts",
    )

    if data.time_left is not None:
        days, hours, minutes, seconds = data.time_left
        time_values = {
            "DAYS": days,
            "HOURS": hours,
            "MINUTES": minutes,
            "SECONDS": seconds,
        }
        for unit, value in time_values.items():
            updated = replace_once(
                updated,
                rf"const\s+int\s+EVENT_DURATION_{unit}\s*=\s*\d+\s*;",
                f"const int EVENT_DURATION_{unit} = {value};",
                f"EVENT_DURATION_{unit}",
            )

    temp_path = cpp_path.with_suffix(cpp_path.suffix + ".tmp")
    temp_path.write_text(updated, encoding="utf-8")
    os.replace(temp_path, cpp_path)


def default_cpp_target(directory: Path) -> Path:
    current_year = datetime.now().year
    candidates = sorted(directory.glob(f"*_{current_year}.cpp"))
    if len(candidates) == 1:
        return candidates[0]
    names = ", ".join(path.name for path in candidates) or "none"
    raise ValueError(
        f"Could not choose a unique *_{current_year}.cpp target (found: {names}). "
        "Pass --cpp explicitly."
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read an ITRTG generator-event screenshot and update the C++ settings."
    )
    parser.add_argument("screenshot", type=Path, help="path to the game screenshot")
    parser.add_argument("--cpp", type=Path, help="C++ optimizer file to update")
    parser.add_argument(
        "--dry-run", action="store_true", help="extract and print values without editing a file"
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="update even when the screenshot's first six item names do not match the target",
    )
    parser.add_argument("--show-ocr", action="store_true", help="print raw OCR text for diagnostics")
    return parser


async def async_main(args: argparse.Namespace) -> int:
    data = await extract_screenshot(args.screenshot)
    print_data(data, show_ocr=args.show_ocr)
    if args.dry_run:
        print("\nDry run: no files changed.")
        return 0

    cpp_path = args.cpp or default_cpp_target(Path.cwd())
    update_cpp(cpp_path, data, force=args.force)
    print(f"\nUpdated {cpp_path.resolve()}")
    if data.time_left is None:
        print("Warning: event time was not recognized, so duration constants were left unchanged.")
    return 0


def main() -> int:
    args = build_parser().parse_args()
    try:
        return asyncio.run(async_main(args))
    except (FileNotFoundError, ImportError, RuntimeError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
