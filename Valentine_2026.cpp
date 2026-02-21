// Original Author & Concept: Syokora
// V2 Update by: Frolf
// added "Busy Times" feature so purchases won't be planned for when you're at work, school, sleeping, etc. Should help make the sims "one and done"
// added multi-thread parallelization to produce multiple solutions and pick the best one
// added additional methods for the optimizer to escape local optima (randomize subsequences, delete upgrades, etc.)
// Made existing methods more thorough to help escape local optima
// Made existing methods more efficient to improve runtime (Should get fairly decent results in ~10 seconds)
// Rehauled variable names, created classes, structs, added modularity, and cleaned up comments to make it easier to maintain for others 
// Improve generalization to make it easier to update for different events;
// Fixed a few non-critical bugs that were causing slower run-times or inadequate searching
// Fixed rare bug where the score can become unbounded and skyrocket
// (WIP) added the option for absolute timestamps for upgrade buy-times rather than timestamps relative to the start of the sim
// (WIP) add parallelization for a isngle solution to produce a near-instant single solution

#include <iostream>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <random> 
#include <array>
#include <unordered_set>

using namespace std;
typedef long long ll;

// USER SETTINGS ---------------------------------------------------------------------
const int EVENT_DURATION_DAYS = 14;
const int EVENT_DURATION_HOURS = 0;
const int EVENT_DURATION_MINUTES = 0;
const int EVENT_DURATION_SECONDS = 0;
const int UNLOCKED_PETS = 56;
const int DLs = 1348;

vector<int> currentLevels = { 
    // Current Production Levels
    0, 1, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0, 
    
    // Current Speed Levels
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    
    0 // Dummy placeholder. Keep 0
};
vector<double> resourceCounts = { 
    // Current resource counts
    0, 500000, 0,    
    0, 0, 0,  
    0/((500.0+DLs)/5.0), 0, 0,
    0, 0, 0/(UNLOCKED_PETS/100.0)                        
};

// Your current upgrade path/the path you want to optimize. If left blank, a random path will be generated.
vector<int> upgradePath = {};
const bool isFullPath = false;
const bool runOptimization = true; // Set false to only see the results and timings of your path
const bool endlessMode = true; //Repeat optimization until *manually* stopped. Prints a path to file if it's better than every other previous path. Only works running locally

// END USER SETTINGS ---------------------------------------------------------------------
// ADVANCED SETTINGS ---------------------------------------------------------------------
const int outputInterval = 5000; //How often improvements are allowed to be printed (milliseconds)

const double EVENT_CURRENCY_WEIGHT = 100e-5; 
const double FREE_EXP_WEIGHT = 6e-5;     
const double PET_STONES_WEIGHT = 3.5e-5;
const double RESEARCH_POINTS_WEIGHT = 2e-5;  
const double GROWTH_WEIGHT = 8e-5;       

const vector<double> busyTimesStart =   {}; // The first hours, relative to the start of the simulation, you'll be unable to purchase upgrades
const vector<double> busyTimesEnd =     {}; // The first hours, relative to the start of the simulation, you'll be able to purchase upgrades again
// END ADVANCED SETTINGS ------------------------------------------------------------------
// PROGRAM SETTINGS -----------------------------------------------------------------------
constexpr array<const char*, 12> resourceNames = {"Love_Bow", "Love_Arrows", "Chocolate", "Cubear", "Rose", "Cake", "FREE_EXP", "Love_Bear", "PET_STONES", "RESEARCH_POINTS", "EVENT_CURRENCY", "GROWTH"};
constexpr double INFINITY_VALUE = (1e100);
constexpr int NUM_RESOURCES = resourceNames.size();
constexpr int NUM_UPGRADES = NUM_RESOURCES * 2;
const int MAX_LEVEL = 71; //Pretend there's a max level for constructing lookup tables. Any number that won't practically be reached is fine to use.
const int MAX_SPEED_LEVEL = 11; // 0-10 is 11 distinct "levels"
constexpr int TOTAL_SECONDS = ((EVENT_DURATION_DAYS)*24*3600+(EVENT_DURATION_HOURS)*3600+(EVENT_DURATION_MINUTES)*60+EVENT_DURATION_SECONDS);
map<int, string> upgradeNames;
array<double, TOTAL_SECONDS> timeNeededSeconds{};

// Rate Cache: [ResourceID][Level][SpeedLevel]
// Cost Cache: [UpgradeID][Current Level]
// We will precompute the production rates for each resource at each level and speed to speed up the simulation.
double RATE_CACHE[NUM_RESOURCES][MAX_LEVEL][MAX_SPEED_LEVEL];

constexpr double cycleTimeMultiplier[12] = {
        1.0 / 3.0,
        1.0,
        1.0 / 3.0,
        1.0 / 5.0,
        1.0 / 5.0,
        1.0 / 5.0,
        1.0 / 2500.0,
        1.0 / 5.0,
        1.0 / 1200.0,
        1.0 / 1200.0,
        1.0 / 5000.0,
        1.0 / 1800.0 
};
    constexpr double speedMultipliers[11] = {
        1.0,
        1.25,
        1.5625,
        1.953125,
        2.44140625,
        3.0517578125,
        3.814697265625,
        4.76837158203125,
        5.960464477539063,
        7.450580596923828,
        9.313225746154785
};
// END PROGRAM SETTINGS -----------------------------------------------------------------
// UTILITY FUNCTIONS --------------------------------------------------------------------
template <typename T>
void printVector(vector<T>& x, ostream& out = cout) {
    for (auto item : x) out << item << ",";
}
class NullBuffer : public std::streambuf {
public:
    int overflow(int c) override { return c; }
};
class NullStream : public std::ostream {
public:
    NullStream() : std::ostream(&m_sb) {}
private:
    NullBuffer m_sb;
};
class Logger {
    int interval;
    mutable chrono::steady_clock::time_point lastLogTime;
    ostream& out;

public:
    Logger(int outputInterval, ostream& output = cout)
        : interval(outputInterval),
          lastLogTime(chrono::steady_clock::now()),
          out(output) {}

    void logImprovement(const string& type, vector<int>& path, const double score) const {
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::milliseconds>(now - lastLogTime).count();

        if (elapsed >= interval) {
            out << "Improved path (" << type << "): \n";
            out << "{";
            printVector(path, out);
            out << "} \n";
            out << "Score: " << score << "\n";
            lastLogTime = now;
        }
    }
};
struct SearchContext {
    Logger logger;
    const vector<double>& resources;
    const vector<int>& levels;
};
struct OptimizationPackage {
    vector<int> path;
    double score;
    mt19937 randomEngine;
    unordered_set<string> deadMoves = {};
};
struct CostEntry {
    int resourceIdx; // The ingredient (0-11)
    double amount;   // The cost
};
vector<CostEntry> COST_CACHE[NUM_UPGRADES+1][MAX_LEVEL];
void nameUpgrades() {
    for (int i = 0; i < NUM_RESOURCES; i++) {
        upgradeNames[i] = string(resourceNames[i]) + "_Level";
        upgradeNames[i + NUM_RESOURCES] = string(resourceNames[i]) + "_Speed";
    }
    upgradeNames[24] = "Complete";
}
vector <int> adjustFullPath(vector<int>& levels){ // This is not thread-safe
    levels[1]--; // Reduce first level because it always starts at 1
    auto it = upgradePath.begin();
    while (it != upgradePath.end()) {
        if (levels[*it] > 0) {
            levels[*it]--;
            it = upgradePath.erase(it);
        } else {
            ++it;
        }
    }
    return levels;
}
void printFormattedResults(vector<int>& path, vector<int>& simulationLevels, vector<double>& simulationResources, double finalScore) {
    cout << "Upgrade Path: \n{";
    printVector(path);
    cout << "}\n";
    cout << "Final Resource Counts: ";
    printVector(simulationResources);
    cout << "\n";
    cout << "Final Upgrade Levels: ";
    printVector(simulationLevels);
    cout << "\n";

    cout << "Event Currency: " << simulationResources[10] << "\n";
    cout << "Free Exp (" << DLs << " DLs): " 
        << simulationResources[6] * (500.0 + DLs) / 5.0 
        << " (" << simulationResources[6] << " levels * cycles)" << "\n";
    cout << "Pet Stones: " << simulationResources[8] << "\n";
    cout << "Research Points: " << simulationResources[9] << "\n";
    cout << "Growth (" << UNLOCKED_PETS << " pets): " 
        << simulationResources[11] * UNLOCKED_PETS / 100.0 
        << " (" << simulationResources[11] << " levels * cycles)" << "\n";
    cout << "Score: " << finalScore << "\n\n";
}
void cacheProductionRates() {
    for (int res = 0; res < NUM_RESOURCES; res++) {
        for (int lvl = 0; lvl < MAX_LEVEL; lvl++) {
            for (int spd = 0; spd < MAX_SPEED_LEVEL; spd++) {
                double rate = (double)lvl * cycleTimeMultiplier[res] * speedMultipliers[spd];
                RATE_CACHE[res][lvl][spd] = rate;
            }
        }
    }
}
vector<double> initializeProductionRates() {
    vector<double> rates(NUM_RESOURCES);
    for (int res = 0; res < NUM_RESOURCES; res++) {
        rates[res] = RATE_CACHE[res][currentLevels[res]][currentLevels[res + NUM_RESOURCES]];
    }
    return rates;
}
void preprocessCosts() {
    for (int type = 0; type < NUM_UPGRADES; type++) {
        int resourceType = type % NUM_RESOURCES;
        if (type == NUM_RESOURCES * 2) resourceType = -1; // Completion upgrade

        for (int lvl = 0; lvl < MAX_LEVEL; lvl++) {
            double nextLvl = (double)(lvl + 1);
            double baseCost = nextLvl * nextLvl * nextLvl * 100.0;
            if (type >= NUM_RESOURCES) baseCost *= 2.5;

            // Define dependencies (Temporary vector to hold multipliers)
            vector<pair<int, double>> dependencies; 
            
            switch (resourceType) {
                case 0:  dependencies = {{1, 2.0}}; break;
                case 1:  dependencies = {{1, 1.0}}; break;
                case 2:  dependencies = {{1, 2.0}}; break;
                case 3:  dependencies = {{0, 1.5}, {1, 1.5}}; break;
                case 4:  dependencies = {{1, 2.0}}; break;
                case 5:  dependencies = {{1, 1.5}, {2, 1.5}}; break;
                case 6:  dependencies = {{3, 5.0}}; break;
                case 7:  dependencies = {{1, 1.5}, {4, 1.5}}; break;
                case 8:  dependencies = {{5, 5.0}}; break;
                case 9:  dependencies = {{3, 2.5}, {4, 2.5}}; break;
                case 10: dependencies = {{3, 3.0}, {4, 3.0}, {5, 3.0}}; break;
                case 11: dependencies = {{5, 2.5}, {7, 2.5}}; break;
                default: break; // Case -1 or others
            }

            // Fill the cache
            for (auto& dep : dependencies) {
                COST_CACHE[type][lvl].push_back(CostEntry{dep.first, baseCost * dep.second});
            }
        }
    }
}

void preprocessBusyTimes(const vector<double>& startHours, const vector<double>& endHours) {
    for (size_t i = 0; i < startHours.size(); ++i) {
        int startSec = static_cast<int>(startHours[i] * 3600.0);
        int endSec = static_cast<int>(endHours[i] * 3600.0);
        startSec = min(startSec, TOTAL_SECONDS - 1);
        endSec = min(endSec, TOTAL_SECONDS - 1);

        for (int s = startSec; s <= endSec; ++s) {
            timeNeededSeconds[s] = endSec - s;
        }
    }
}
inline double additionalTimeNeeded(double expectedTimeSeconds) {
    int idx = static_cast<int>(expectedTimeSeconds);
    if (idx >= TOTAL_SECONDS) return 0.0;
    return timeNeededSeconds[idx];
}
vector <int> generateRandomPath(int length = TOTAL_SECONDS/3600) {
    vector<int> randomPath = {};
    for (int i = 0; i < length; i++) {
        randomPath.push_back(rand() % NUM_RESOURCES + NUM_RESOURCES * (rand() % 2));
    }
    randomPath.push_back(NUM_RESOURCES * 2);
    return randomPath;
}
void readoutUpgrade(int upgradeType, vector<int>& levels, int elapsedSeconds) {
    cout << upgradeNames[upgradeType] << " " << levels[upgradeType] 
        << " " << (int)elapsedSeconds/24/3600 << " days, " 
        << (int)elapsedSeconds/3600%24 << " hours, " 
        << (int)elapsedSeconds/60%60 << " minutes" << "\n";
}
// END UTILITY FUNCTIONS ----------------------------------------------------------------
// ALGORITHM FUNCTIONS ------------------------------------------------------------------
double performUpgrade(vector<int>& levels, vector<double>& resources, vector<double>& currentRates, int upgradeType, double& remainingTime) {
    
    int currentLevel = levels[upgradeType];
    
    const vector<CostEntry>& costs = COST_CACHE[upgradeType][currentLevel];
    double timeNeeded = 0;
    double timeElapsed = TOTAL_SECONDS - remainingTime;
    for (const auto& entry : costs) {
        double neededResource = entry.amount - resources[entry.resourceIdx];
        double rate = RATE_CACHE[entry.resourceIdx][levels[entry.resourceIdx]][levels[entry.resourceIdx + NUM_RESOURCES]];
        if (neededResource <= 0) continue;
        if (rate == 0) {
            timeNeeded = INFINITY_VALUE;
            break;
        }
        timeNeeded = max(timeNeeded, neededResource / rate);
    }

    int busyLookupIndex = static_cast<int>(timeElapsed + timeNeeded);
    if (0 <= busyLookupIndex && busyLookupIndex < TOTAL_SECONDS){
        timeNeeded += timeNeededSeconds[busyLookupIndex];
    };

    if (timeNeeded >= remainingTime || upgradeType == (2 * NUM_RESOURCES)) {
        timeNeeded = remainingTime;
    }

    resources[0] += currentRates[0] * timeNeeded;
    resources[1] += currentRates[1] * timeNeeded;
    resources[2] += currentRates[2] * timeNeeded;
    resources[3] += currentRates[3] * timeNeeded;
    resources[4] += currentRates[4] * timeNeeded;
    resources[5] += currentRates[5] * timeNeeded;
    resources[6] += currentRates[6] * timeNeeded;
    resources[7] += currentRates[7] * timeNeeded;
    resources[8] += currentRates[8] * timeNeeded;
    resources[9] += currentRates[9] * timeNeeded;
    resources[10] += currentRates[10] * timeNeeded;
    resources[11] += currentRates[11] * timeNeeded;

    if (timeNeeded < remainingTime) {
        for (const auto & entry : costs) {
            resources[entry.resourceIdx] -= entry.amount;
        }
        levels[upgradeType]++;
        int resChanged = upgradeType % NUM_RESOURCES;
        int lvl = levels[resChanged];
        int spd = levels[resChanged + NUM_RESOURCES];
        currentRates[resChanged] = RATE_CACHE[resChanged][lvl][spd];
    }
    
    return timeNeeded;
}
double simulateUpgradePath(vector<int>& path, vector<int>& levels, vector<double>& resources, bool display = false) {
    double time = TOTAL_SECONDS;
    thread_local vector<double> currentRates;
    currentRates = initializeProductionRates();

    for (auto upgradeType : path) {
        if (time < 1e-3) return 0;
        if (upgradeType >= NUM_RESOURCES && levels[upgradeType] >= 10) return INFINITY_VALUE;
        double timeTaken = performUpgrade(levels, resources, currentRates, upgradeType, time);
        time -= timeTaken;
        
        if (display) {
            int elapsedSeconds = TOTAL_SECONDS - time;
            readoutUpgrade(upgradeType, levels, elapsedSeconds);
        }
    }
    return time;
}
double calculateScore(vector<double>& resources, bool display = false) {
    double score = 0;

    for (int i = 0; i < NUM_RESOURCES; i++) {
        score += resources[i] * 1e-15;
    }

    score += (min(resources[10], 10000.0) + max(0.0, (resources[10] - 10000)) * 0.01) * (EVENT_CURRENCY_WEIGHT);
    score += resources[6] * (FREE_EXP_WEIGHT);          // Free EXP
    score += resources[8] * (PET_STONES_WEIGHT);        // Pet Stones
    score += resources[9] * (RESEARCH_POINTS_WEIGHT);   // Research Points
    score += resources[11] * (GROWTH_WEIGHT);           // Growth
    
    return score;
}
double evaluatePath(vector<int>& path, const SearchContext& context){
    thread_local vector<double> testResources = context.resources;
    thread_local vector<int> testLevels = context.levels;
    testResources = context.resources;
    testLevels = context.levels;
    simulateUpgradePath(path, testLevels, testResources);
    return calculateScore(testResources);
}
double calculateFinalPath(vector<int>& path){
    vector<int>     simulationLevels(currentLevels);
    vector<double>  simulationResources(resourceCounts);

    simulateUpgradePath(path, simulationLevels, simulationResources, true);
    double simulationScore = calculateScore(simulationResources);
    printFormattedResults(path, simulationLevels, simulationResources, simulationScore);
    return simulationScore;
}
bool tryInsertUpgrade(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size();

    uniform_int_distribution<> positionDist(0, pathLength);
    int startPosition = positionDist(package.randomEngine);
    const int maxTypes = (NUM_RESOURCES * 2);

    for (int i = 0; i < pathLength; i++) {
        int modulatedInsertPosition = (i + startPosition) % pathLength;
        package.path.insert(package.path.begin() + modulatedInsertPosition, 0);

        int startingUpgradeType = package.randomEngine() % (NUM_RESOURCES * 2);
        for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {

            int modulatedUpgradeType = (upgradeType + startingUpgradeType) % maxTypes;
            package.path[modulatedInsertPosition] = modulatedUpgradeType;

            double testScore = evaluatePath(package.path, context);

            if (testScore > package.score) {
                package.score = testScore;
                context.logger.logImprovement("Insert", package.path, package.score);
                return true;
            }
        }
        package.path.erase(package.path.begin() + modulatedInsertPosition);
    }
    package.deadMoves.insert("Insert");
    return false;
}
bool tryRemoveUpgrade(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size() - 1;

    uniform_int_distribution<> swapDist(0, pathLength - 2);
    int startPos = swapDist(package.randomEngine);

    for (int i = 0; i < pathLength; i++) {
        int removePos = (i + startPos) % (pathLength);
        int removedUpgrade = package.path[removePos];

        package.path.erase(package.path.begin() + removePos);
        double testScore = evaluatePath(package.path, context);

        if (testScore >= package.score) {
            package.score = testScore;
            context.logger.logImprovement("Remove", package.path, package.score);
            return true;
        }
        package.path.insert(package.path.begin() + removePos, removedUpgrade);
    }
    package.deadMoves.insert("Remove");
    return false;

}
bool tryReplaceUpgrade(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size();

    uniform_int_distribution<> positionDist(0, pathLength - 1); // -1 because we access index, not insert
    int startPosition = positionDist(package.randomEngine);
    const int maxTypes = (NUM_RESOURCES * 2);

    for (int i = 0; i < pathLength; i++) {

        int targetIndex = (i + startPosition) % pathLength;
        int originalType = package.path[targetIndex];

        int startingUpgradeType = package.randomEngine() % (NUM_RESOURCES * 2);
        
        for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {
            int modulatedUpgradeType = (upgradeType + startingUpgradeType) % maxTypes;

            if (modulatedUpgradeType == originalType) continue;

            package.path[targetIndex] = modulatedUpgradeType;

            double testScore = evaluatePath(package.path, context);

            if (testScore > package.score) {
                package.score = testScore;
                context.logger.logImprovement("Replace", package.path, package.score);
                return true;
            }
        }
        package.path[targetIndex] = originalType;
    }
    package.deadMoves.insert("Replace");
    return false;
}
bool trySwapUpgrades(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size() - 1;
    double testScore;

    uniform_int_distribution<> swapDist(0, pathLength - 2);
    int startPos = swapDist(package.randomEngine);

    for (int i2 = 0; i2 < pathLength - 1; i2++) { 
        for (int j2 = i2 + 1; j2 < pathLength - 1; j2++) { 
            int i = (i2 + startPos) % (pathLength - 1);
            int j = (j2 + startPos) % (pathLength - 1);

            if (package.path[i] == package.path[j]) continue;
            swap(package.path[i], package.path[j]);

            testScore = evaluatePath(package.path, context);

            if (testScore > package.score) {
                package.score = testScore;
                context.logger.logImprovement("Swap", package.path, package.score);
                return true;
            }
            swap(package.path[i], package.path[j]);
        }
    }
    package.deadMoves.insert("Swap");
    return false;
}
bool tryRotateSubsequences(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size() - 1;
    int maxIndex = pathLength - 1;
    double testScore;

    // Pick a Random Subsequence Rage [rangeStart, rangeEnd]
    // We need at least 3 items to perform meaningful rotation (2 items is just a swap)
    uniform_int_distribution<> startDist(0, maxIndex - 2);
    int rangeStart = startDist(package.randomEngine);

    uniform_int_distribution<> endDist(rangeStart + 2, maxIndex);
    int rangeEnd = endDist(package.randomEngine);

    // Iterators for the range [first, last)
    // Note: rangeEnd is inclusive index, so iterator is +1
    auto rangeBeginIt = package.path.begin() + rangeStart;
    auto rangeEndIt = package.path.begin() + rangeEnd + 1;
    
    int subSeqLength = rangeEnd - rangeStart;
    for (int attempt = 0; attempt < subSeqLength; attempt++) {

        // Calculate shift logic to try small rotations first (Left 1, Right 1, Left 2, Right 2...)
        int shiftAmount = (attempt + 2) / 2;
        bool rotateLeft = (attempt % 2 == 0);

        if(rotateLeft)  rotate(rangeBeginIt, rangeBeginIt + shiftAmount, rangeEndIt); // Left rotation
        else            rotate(rangeBeginIt, rangeEndIt - shiftAmount, rangeEndIt); // Right rotation

        testScore = evaluatePath(package.path, context);
        if (testScore > package.score) {
            package.score = testScore;
            context.logger.logImprovement("Rotation", package.path, package.score);
            return true;
        }
        //To undo Left(X), we Rotate Right(X)
        //To undo Right(X), we Rotate Left(X)
        if(rotateLeft)  rotate(rangeBeginIt, rangeEndIt - shiftAmount, rangeEndIt);
        else            rotate(rangeBeginIt, rangeBeginIt + shiftAmount, rangeEndIt);
    }
    return false;
}
bool exhaustRotateSubsequences(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size() - 1;
    int maxIndex = pathLength - 1;
    double testScore;

    uniform_int_distribution<> startOffsetDist(0, maxIndex - 2);
    int startOffset = startOffsetDist(package.randomEngine);

    for (int i = 0; i < maxIndex - 1; i++){
        int rangeStart = (startOffset + i) % (maxIndex - 1);

        int remainingLength = maxIndex - (rangeStart + 2);
        uniform_int_distribution<> endOffsetDist(0, max(0, remainingLength)); 
        int endOffset = endOffsetDist(package.randomEngine);
        for (int j = 0; j < remainingLength; j++){ 
            int offsetJ = (endOffset + j) % (remainingLength + 1);
            int rangeEnd = rangeStart + 2 + offsetJ;

            auto rangeBeginIt = package.path.begin() + rangeStart;
            auto rangeEndIt = package.path.begin() + rangeEnd + 1;
            int subSeqLength = rangeEnd - rangeStart;

            for (int attempt = 0; attempt < subSeqLength; attempt++) {

                // Calculate shift logic to try small rotations first (Left 1, Right 1, Left 2, Right 2...)
                int shiftAmount = (attempt + 2) / 2;
                bool rotateLeft = (attempt % 2 == 0);

                if(rotateLeft)  rotate(rangeBeginIt, rangeBeginIt + shiftAmount, rangeEndIt); // Left rotation
                else            rotate(rangeBeginIt, rangeEndIt - shiftAmount, rangeEndIt); // Right rotation

                double testScore = evaluatePath(package.path, context);

                if (testScore > package.score) {
                    package.score = testScore;
                    context.logger.logImprovement("Rotation", package.path, package.score);
                    return true;
                }
                //To undo Left(X), we Rotate Right(X)
                //To undo Right(X), we Rotate Left(X)
                if(rotateLeft)  rotate(rangeBeginIt, rangeEndIt - shiftAmount, rangeEndIt);
                else            rotate(rangeBeginIt, rangeBeginIt + shiftAmount, rangeEndIt);
            }
        }
    }
    package.deadMoves.insert("Rotate");
    return false;
}
void optimizeUpgradePath(OptimizationPackage& package, SearchContext& context, const int maxIterations = 10000) {

    random_device seed;
    mt19937 randomEngine(seed());

    int iterationCount = 0;
    int noImprovementStreak = 0;
    int InsertCount = 0, RemoveCount = 0, SwapCount = 0, RotationCount = 0, ReplaceCount = 0;
    double InsertRatio = 0, RemoveRatio = 0, SwapRatio = 0, RotationRatio = 0, ReplaceRatio = 0;
    while (noImprovementStreak < maxIterations) {
        iterationCount++;
        bool improved = false;

        int strategy = iterationCount % 100;

        if(package.deadMoves.count("Rotate")){
            break;
        }
        else if (package.deadMoves.count("Insert") && package.deadMoves.count("Remove") && package.deadMoves.count("Replace") && package.deadMoves.count("Swap"))
                                                                         {improved = tryRotateSubsequences(package, context); if (improved) RotationCount++;}
        else if (strategy < 15 && !package.deadMoves.count("Insert"))    {improved = tryInsertUpgrade(package, context); if (improved) InsertCount++;}
        else if (strategy < 30 && !package.deadMoves.count("Replace"))   {improved = tryReplaceUpgrade(package, context); if (improved) ReplaceCount++;}
        else if (strategy < 45 && !package.deadMoves.count("Remove"))    {improved = tryRemoveUpgrade(package, context); if (improved) RemoveCount++;}
        else if (strategy < 47)                                          {improved = tryRotateSubsequences(package, context); if (improved) RotationCount++;}
        else if (!package.deadMoves.count("Swap"))                       {improved = trySwapUpgrades(package, context); if (improved) SwapCount++;}
        
        if (improved) {
            noImprovementStreak = 0;
            package.deadMoves.clear();
            continue;
        } 
        else noImprovementStreak++;
    }
    //Data collection for analysis. Realistically needs its own manager.
    InsertRatio = (double)InsertCount / iterationCount;
    ReplaceRatio = (double)ReplaceCount / iterationCount;
    RemoveRatio = (double)RemoveCount / iterationCount;
    SwapRatio = (double)SwapCount / iterationCount;
    RotationRatio = (double)RotationCount / iterationCount;
    cout << "Insert Count: " << InsertCount << "\n";
    cout << "Replace Count: " << ReplaceCount << "\n";
    cout << "Remove Count: " << RemoveCount << "\n";
    cout << "Swap Count: " << SwapCount << "\n";
    cout << "Rotation Count: " << RotationCount << "\n";
    cout << "Insert ratio: " << InsertRatio << "\n";
    cout << "Replace ratio: " << ReplaceRatio << "\n";
    cout << "Remove ratio: " << RemoveRatio << "\n";
    cout << "Swap ratio: " << SwapRatio << "\n";
    cout << "Rotation ratio: " << RotationRatio << "\n";
}

int main() {
    upgradePath.reserve(500);
    vector<int> levelsCopy(currentLevels);
    double finalScore;
    double UltraScore = 0;
    nameUpgrades();
    preprocessBusyTimes(busyTimesStart, busyTimesEnd);
    cacheProductionRates();
    preprocessCosts();
    ofstream MyFile;
    MyFile.open("OptimizedPaths.txt", ios::app);
    while (endlessMode) {
        if (upgradePath.empty()) {          
            upgradePath = generateRandomPath();
        }

        if (isFullPath) {
            levelsCopy = adjustFullPath(levelsCopy);
        }

        double initial_score = calculateFinalPath(upgradePath);
        random_device seed;
        mt19937 randomEngine(seed());
        Logger logger(outputInterval);
        SearchContext context{logger, resourceCounts, currentLevels};
        OptimizationPackage package = {upgradePath, initial_score, move(randomEngine)};
        optimizeUpgradePath(package, context);
        upgradePath = move(package.path);
        
        finalScore =calculateFinalPath(upgradePath);
        if (finalScore >= UltraScore) {
            UltraScore = finalScore;
            printVector(upgradePath, MyFile);
            MyFile << endl;
            MyFile << "Final score: " << finalScore << endl << endl;
        }
        upgradePath = {};
    }

    if (upgradePath.empty()) {          
        upgradePath = generateRandomPath();
    }
    if (isFullPath) {
        levelsCopy = adjustFullPath(levelsCopy);
    }

    double initial_score = calculateFinalPath(upgradePath);

    if (runOptimization) {
        random_device seed;
        mt19937 randomEngine(seed());
        Logger logger(outputInterval);
        SearchContext context{logger, resourceCounts, currentLevels};
        OptimizationPackage package = {upgradePath, initial_score, move(randomEngine)};
        optimizeUpgradePath(package, context);
        upgradePath = move(package.path);
    }
    finalScore = calculateFinalPath(upgradePath);
    
    return 0;
}   

