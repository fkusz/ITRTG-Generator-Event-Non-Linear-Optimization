# ITRTG Generator Event Non-Linear Optimization

A highly optimized C++ tool to calculate and simulate the optimal upgrade path for **Idling to Rule the Gods (ITRTG)** generator events. 

This script uses non-linear optimization techniques to find the best upgrade strategy based on your game variables, saving you time and maximizing your event rewards.

---

## Quick Start (Online Compiler)

Don't want to mess with local environments? You can run the optimizer directly in your browser!

1. **Open an Online Compiler**: Use [C++ Shell](https://cpp.sh/) (or any other C++ compiler supporting **C++23** or newer).
   (Note: some online compilers, such as [Programiz C++ Compiler](https://www.programiz.com/cpp-programming/online-compiler/), do not allow you to choose your compiler version. These will not work for this code.)
3. **Copy the Script**: Paste the code from `EVENT_NAME_YEAR.cpp` into the compiler.
4. **Update `#USER SETTINGS`**:
   - `unlockedPets`: Number of unlocked pets.
   - `DLs`: Number of Dungeon Levels per your stats page.
   - `timeUntilEnd`: Time remaining until the event ends.
   *(Make sure to update these as needed throughout the event!)*
5. **Run the Code**: The first run will time itself out or stop when the score no longer changes drastically. Each block of output shows your **upgrade path** (a long list of comma-separated values) and the **projected score** underneath.
6. **Resume Optimization**: If you want to optimize further, copy the comma-separated values into the `upgradePath` variable (around line 59) in the script. The sim will resume from where it left off.
7. **Get Your Path**: The sim will naturally finish after 10,000 consecutive failed improvements and print a plain, readable upgrade plan. 
   - *Skipping the wait:* If you are happy with a path early, paste the CSV back into `upgradePath` (line 59), set the appropriate flag to `false` (line 62), and run it again to instantly get a readable plan.

> **Note:** The time values printed for each step are relative to the *start* of the script. If you see `28min` followed by `38min`, you wait 10 minutes between purchasing those upgrades.

---

## Running Locally

For maximum performance or endless overnight optimization, run the script locally.

1. **Prerequisites**: Ensure you have an up-to-date compiler supporting **C++23**. Older standards may throw errors. You can use standard tools like `g++`, `clang`, MSVC, or build using the provided `CMakeLists.txt`.
2. **Endless Mode**: Toggle the endless option (around line 63) if you want the program to run for hours finding the absolute best possible path. Otherwise, it ends after 10,000 consecutive failures to improve.
3. **Retrieving Results**: If you manually stop the sim before it finishes naturally, take the best path printed, paste it into `upgradePath` (line 59), recompile, and run again. It will instantly print the readable upgrade guide.

---

## Advanced Settings

Configure these settings inside the script to fine-tune the sim's behavior:

- **Output Interval**: Printing results to the console is slow. Adjust the output interval to reduce console clutter and speed up the sim's execution. Recommended: `3000ms - 5000ms` (3-5 seconds). Very small values will make the program sluggish.
- **Weights**: Prioritize certain event resources over others. Change the weights to score specific outputs more heavily. 
  - *Keep event currency decently higher than other weights, or the sim might ignore it.* Higher relative weights force the sim to produce more of that specific resource.

---

## Scheduling "Busy Times"

The optimizer can account for work meetings, sleeping, or anytime you are away from the game! Pre-schedule your unavailability so the sim can optimize *around* your schedule.

Like output times, these arrays are **relative to the start of your simulation** (in hours).

### Example: Sleeping
If you are starting the simulation at exactly **12:00 PM (Noon)** and plan to sleep from **10:00 PM to 6:00 AM (22:00 - 06:00)**:

```cpp
busyTimesStart = {10, 34, 58, 82};
busyTimesEnd   = {18, 42, 66, 90};
```
* **Explanation:** You are telling the optimizer you'll step away starting in `10` hours (10:00 PM), and you'll be back in `18` hours (6:00 AM). Simply add `24` to each value to repeat it for subsequent days.
