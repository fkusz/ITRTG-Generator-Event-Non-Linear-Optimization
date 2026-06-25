import asyncio
import re
import os
import sys

try:
    from winsdk.windows.media.ocr import OcrEngine
    from winsdk.windows.graphics.imaging import BitmapDecoder
    from winsdk.windows.storage import StorageFile, FileAccessMode
except ImportError:
    print("Error: Missing required Python packages.")
    print("Please install them by running: pip install winsdk pillow")
    sys.exit(1)

async def extract_text(file_path):
    file_path = os.path.abspath(file_path)
    if not os.path.exists(file_path):
        print(f"Error: Could not find image at {file_path}")
        sys.exit(1)
        
    file = await StorageFile.get_file_from_path_async(file_path)
    stream = await file.open_async(FileAccessMode.READ)
    decoder = await BitmapDecoder.create_async(stream)
    software_bitmap = await decoder.get_software_bitmap_async()
    
    engine = OcrEngine.try_create_from_user_profile_languages()
    if not engine:
        print("Error: Could not obtain OCR Engine for your user profile language.")
        return ""
        
    result = await engine.recognize_async(software_bitmap)
    return result.text

def process_and_update(text):
    # The Windows OCR usually reads left-to-right, row-by-row.
    # Therefore, we just dynamically find all Levels, Speeds, and Counts.
    levels = [int(x) for x in re.findall(r'(?:Level|Leve1|Lvl).*?(?::|\s)\s*(\d+)', text, re.IGNORECASE)]
    speeds = [int(x) for x in re.findall(r'(?:Speed|Spd).*?(?::|\s)\s*(\d+)', text, re.IGNORECASE)]
    
    # Counts can contain commas
    counts_raw = re.findall(r'(?:Count|Cunt).*?(?::|\s)\s*([\d,]+)', text, re.IGNORECASE)
    counts = [int(x.replace(',', '')) for x in counts_raw]

    if len(levels) != 12 or len(speeds) != 12 or len(counts) != 12:
        print("Warning: OCR did not find exactly 12 of each stat.")
        print(f"Found {len(levels)} Levels, {len(speeds)} Speeds, {len(counts)} Counts.")
        print("Ensure the screenshot is clear and captures all 12 items.")
        print("Continuing with best effort...\n")
        
        # Pad with 0s if missing
        while len(levels) < 12: levels.append(0)
        while len(speeds) < 12: speeds.append(0)
        while len(counts) < 12: counts.append(0)

    print("--- Extracted Values ---")
    print(f"Levels: {levels[:12]}")
    print(f"Speeds: {speeds[:12]}")
    print(f"Counts: {counts[:12]}")
    print("------------------------\n")

    cpp_file = "Easter_2026.cpp"
    if not os.path.exists(cpp_file):
        print(f"Error: Could not find {cpp_file} in current directory.")
        return

    with open(cpp_file, 'r', encoding='utf-8') as f:
        content = f.read()

    # Format the replacement arrays
    levels_str = f"""array<int, 25> currentLevels = {{ 
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

    # Format resource counts (re-applying proper mathematical expansions for indices 6, 9, 11)
    counts_str = f"""array<double, 12> resourceCounts = {{ 
    // Current resource counts
    {counts[0]}, {counts[1]}, {counts[2]}, 
    {counts[3]}, {counts[4]}, {counts[5]}, 
    {counts[6]}/((500.0+DLs)/5.0), {counts[7]}, {counts[8]},
    {counts[9]}/(0.5+ML/100.0), {counts[10]}, {counts[11]}/(UNLOCKED_PETS/100.0)
}};"""

    # Do Regex replace
    content = re.sub(r'array<int, 25> currentLevels.*?};', levels_str, content, flags=re.DOTALL)
    content = re.sub(r'array<double, 12> resourceCounts.*?};', counts_str, content, flags=re.DOTALL)

    with open(cpp_file, 'w', encoding='utf-8') as f:
        f.write(content)
        
    print(f"Successfully updated {cpp_file}!")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python update_from_screenshot.py <screenshot_path>")
        sys.exit(1)
        
    image_path = sys.argv[1]
    
    # Run async function using standard asyncio
    text_data = asyncio.run(extract_text(image_path))
    
    if text_data:
        process_and_update(text_data)
