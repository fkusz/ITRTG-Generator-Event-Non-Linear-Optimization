#ifndef ENABLE_PARALLEL_SWAP
#define ENABLE_PARALLEL_SWAP 1 // Set to 0 if your compiler has trouble linking std::thread.
#endif

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
#if ENABLE_PARALLEL_SWAP
#include <thread>
#include <atomic>
#endif

using namespace std;
typedef long long ll;

// USER SETTINGS ---------------------------------------------------------------------
const int EVENT_DURATION_DAYS = 14;
const int EVENT_DURATION_HOURS = 0;
const int EVENT_DURATION_MINUTES = 0;
const int EVENT_DURATION_SECONDS = 0;
const int UNLOCKED_PETS = 47;
const int DLs = 3070;
const int AL = 211;

array<int, 25> currentLevels = { 
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
array<double, 12> resourceCounts = { 
    // Current resource counts
    0, 500000, 0, 
    0, 0, 0, 
    0/((500.0+DLs)/5.0), 0, 0,
    0/(0.5+AL/100.0), 0, 0/(UNLOCKED_PETS/100.0)
};

// Your current upgrade path/the path you want to optimize. If left blank, a random path will be generated.
vector<int> upgradePath = {};
const bool isFullPath = false;
const bool runOptimization = true; // Set false to only see the results and timings of your path
const bool endlessMode = false; //Repeat optimization until *manually* stopped. Prints a path to file if it's better than every other previous path. Only works running locally

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
constexpr array<const char*, 12> resourceNames = {"Blue", "White", "Mixed", "Blue/White", "Red", "Striped", "FREE_EXP", "Red/White", "PET_STONES", "RESEARCH_POINTS", "EVENT_CURRENCY", "GROWTH"};
constexpr double INFINITY_VALUE = (1e100);
constexpr int NUM_RESOURCES = resourceNames.size();
constexpr int NUM_UPGRADES = NUM_RESOURCES * 2;
const int MAX_LEVEL = 71; //Pretend there's a max level for constructing lookup tables. Any number that won't practically be reached is fine to use.
const int MAX_SPEED_LEVEL = 11; // 0-10 is 11 distinct "levels"
#if ENABLE_PARALLEL_SWAP
constexpr long long PARALLEL_SWAP_MIN_CANDIDATES = 20000;
constexpr int PARALLEL_SWAP_MAX_THREADS = 8;
#endif
constexpr int TOTAL_SECONDS = ((EVENT_DURATION_DAYS)*24*3600+(EVENT_DURATION_HOURS)*3600+(EVENT_DURATION_MINUTES)*60+EVENT_DURATION_SECONDS);
map<int, string> upgradeNames;
array<double, TOTAL_SECONDS + 1> timeNeededSeconds{};

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
void printVector(const auto& x, ostream& out = cout) {
    for (const auto& item : x) out << item << ",";
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
    const array<double, NUM_RESOURCES>& resources;
    const array<int, NUM_UPGRADES + 1>& levels;
};
struct SimSnapshot {
    array<double, NUM_RESOURCES> resources;
    array<double, NUM_RESOURCES> rates;
    array<int, NUM_UPGRADES + 1> levels;
    double remainingTime;
};
struct OptimizationPackage {
    vector<int> path;
    double score;
    mt19937 randomEngine;
    vector<SimSnapshot> snapshots;
    bool snapshotsValid = false;
    bool deadInsert = false;
    bool deadRemove = false;
    bool deadReplace = false;
    bool deadSwap = false;
    bool deadRotate = false;
};
struct CostEntry {
    int resourceIdx; // The ingredient (0-11)
    double amount;   // The cost
};
struct CostList {
    array<CostEntry, 3> entries;
    int count = 0;
};

CostList COST_CACHE[NUM_UPGRADES+1][MAX_LEVEL];
void nameUpgrades() {
    for (int i = 0; i < NUM_RESOURCES; i++) {
        upgradeNames[i] = string(resourceNames[i]) + "_Level";
        upgradeNames[i + NUM_RESOURCES] = string(resourceNames[i]) + "_Speed";
    }
    upgradeNames[24] = "Complete";
}
array<int,NUM_UPGRADES + 1> adjustFullPath(array<int,NUM_UPGRADES + 1>& levels){ // This is not thread-safe
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
void printFormattedResults(vector<int>& path, array<int, NUM_UPGRADES + 1>& simulationLevels, array<double, NUM_RESOURCES>& simulationResources, double finalScore) {
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
    cout << "Research Points: " << simulationResources[9] * (0.5 + AL/100.0) << "\n";
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
array<double, NUM_RESOURCES> initializeProductionRates(const array<int, NUM_UPGRADES + 1>& levels = currentLevels) {
    array<double, NUM_RESOURCES> rates;
    for (int res = 0; res < NUM_RESOURCES; res++) {
        rates[res] = RATE_CACHE[res][levels[res]][levels[res + NUM_RESOURCES]];
    }
    return rates;
}
void preprocessCosts() {
    for (int type = 0; type <= NUM_UPGRADES; type++) { 
        int resourceType = type % NUM_RESOURCES;
        if (type == NUM_RESOURCES * 2) resourceType = -1; // Completion upgrade

        for (int lvl = 0; lvl < MAX_LEVEL; lvl++) {
            
            CostList& costList = COST_CACHE[type][lvl];
            costList.count = 0;

            // Give the completion upgrade an impossible cost to ensure it ends performUpgrade()
            if (resourceType == -1) {
                costList.entries[0] = CostEntry{1, INFINITY_VALUE};
                costList.count = 1;
                continue;
            }

            double nextLvl = (double)(lvl + 1);
            double baseCost = nextLvl * nextLvl * nextLvl * 100.0;
            if (type >= NUM_RESOURCES) baseCost *= 2.5;

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
            }

            for (size_t i = 0; i < dependencies.size(); ++i) {
                costList.entries[i] = CostEntry{dependencies[i].first, baseCost * dependencies[i].second};
            }
            costList.count = static_cast<int>(dependencies.size());
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
vector <int> generateRandomPath(int length = TOTAL_SECONDS/3600) {
    vector<int> randomPath = {};
    for (int i = 0; i < length; i++) {
        randomPath.push_back(rand() % NUM_RESOURCES + NUM_RESOURCES * (rand() % 2));
    }
    randomPath.push_back(NUM_RESOURCES * 2);
    return randomPath;
}
void readoutUpgrade(int upgradeType, array<int,NUM_UPGRADES + 1>& levels, int elapsedSeconds) {
    cout << upgradeNames[upgradeType] << " " << levels[upgradeType] 
        << " " << (int)elapsedSeconds/24/3600 << " days, " 
        << (int)elapsedSeconds/3600%24 << " hours, " 
        << (int)elapsedSeconds/60%60 << " minutes, "
        << (int)elapsedSeconds%60 << " seconds" << "\n";
}
// END UTILITY FUNCTIONS ----------------------------------------------------------------
// ALGORITHM FUNCTIONS ------------------------------------------------------------------
double performUpgrade(array<int, NUM_UPGRADES + 1>& levels, 
                      array<double, NUM_RESOURCES>& resources, 
                      array<double, NUM_RESOURCES>& currentRates, 
                      int upgradeType, 
                      double& remainingTime) {
    
    int currentLevel = levels[upgradeType];
    const CostList& costList = COST_CACHE[upgradeType][currentLevel];
    
    double timeNeeded = 0;
    double timeElapsed = TOTAL_SECONDS - remainingTime;
    
    // Flattened cache loop - Zero heap indirection!
    for (int i = 0; i < costList.count; ++i) {
        const CostEntry& entry = costList.entries[i];
        double neededResource = entry.amount - resources[entry.resourceIdx];
        if (neededResource <= 0) continue;
        double rate = currentRates[entry.resourceIdx];
        if (rate == 0) {
            timeNeeded = remainingTime;
            break;
        }
        timeNeeded = max(timeNeeded, neededResource / rate);
    }

    timeNeeded = min(timeNeeded, remainingTime);
    int busyLookupIndex = static_cast<int>(timeElapsed + timeNeeded);
    timeNeeded += timeNeededSeconds[busyLookupIndex];

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
        for (int i = 0; i < costList.count; ++i) {
            const CostEntry& entry = costList.entries[i];
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
SimSnapshot makeInitialSnapshot(const SearchContext& context) {
    SimSnapshot snapshot;
    snapshot.resources = context.resources;
    snapshot.levels = context.levels;
    snapshot.rates = initializeProductionRates(snapshot.levels);
    snapshot.remainingTime = TOTAL_SECONDS;
    return snapshot;
}
double simulateUpgradePathFrom(const vector<int>& path, SimSnapshot& snapshot, int startIndex = 0, bool display = false) {
    for (int pathIndex = startIndex; pathIndex < static_cast<int>(path.size()); pathIndex++) {
        if (snapshot.remainingTime < 1e-3) return 0;

        int upgradeType = path[pathIndex];
        if (upgradeType >= NUM_RESOURCES && snapshot.levels[upgradeType] >= 10) {
            snapshot.remainingTime = 0;
            return INFINITY_VALUE;
        }

        double timeTaken = performUpgrade(snapshot.levels, snapshot.resources, snapshot.rates, upgradeType, snapshot.remainingTime);
        snapshot.remainingTime -= timeTaken;
        
        if (display) {
            int elapsedSeconds = TOTAL_SECONDS - snapshot.remainingTime;
            readoutUpgrade(upgradeType, snapshot.levels, elapsedSeconds);
        }
    }
    return snapshot.remainingTime;
}
double simulateUpgradePath(const vector<int>& path, array<int,NUM_UPGRADES + 1>& levels, array<double, NUM_RESOURCES>& resources, bool display = false) {
    SimSnapshot snapshot;
    snapshot.resources = resources;
    snapshot.levels = levels;
    snapshot.rates = initializeProductionRates(snapshot.levels);
    snapshot.remainingTime = TOTAL_SECONDS;

    double remainingTime = simulateUpgradePathFrom(path, snapshot, 0, display);
    resources = snapshot.resources;
    levels = snapshot.levels;
    return remainingTime;
}
double calculateScore(array<double, NUM_RESOURCES>& resources, bool display = false) {
    double score = 0;

    for (int i = 0; i < NUM_RESOURCES; i++) {
        score += resources[i] * 1e-15; // Give a very small score for each resource to incentivize early investment. Otherwise sim breaks for event start
    }

    score += (min(resources[10], 10000.0) + max(0.0, (resources[10] - 10000)) * 0.005) * (EVENT_CURRENCY_WEIGHT);
    score += resources[6] * (FREE_EXP_WEIGHT);
    score += resources[8] * (PET_STONES_WEIGHT);
    score += resources[9] * (RESEARCH_POINTS_WEIGHT);
    score += resources[11] * (GROWTH_WEIGHT);
    
    return score;
}
void buildSnapshots(const vector<int>& path, const SearchContext& context, vector<SimSnapshot>& snapshots) {
    snapshots.clear();
    snapshots.reserve(path.size() + 1);

    SimSnapshot snapshot = makeInitialSnapshot(context);
    snapshots.push_back(snapshot);

    for (int pathIndex = 0; pathIndex < static_cast<int>(path.size()); pathIndex++) {
        if (snapshot.remainingTime >= 1e-3) {
            int upgradeType = path[pathIndex];
            if (upgradeType >= NUM_RESOURCES && snapshot.levels[upgradeType] >= 10) {
                snapshot.remainingTime = 0;
            } else {
                double timeTaken = performUpgrade(snapshot.levels, snapshot.resources, snapshot.rates, upgradeType, snapshot.remainingTime);
                snapshot.remainingTime -= timeTaken;
            }
        }
        snapshots.push_back(snapshot);
    }
}
void ensureSnapshots(OptimizationPackage& package, const SearchContext& context) {
    if (!package.snapshotsValid) {
        buildSnapshots(package.path, context, package.snapshots);
        package.snapshotsValid = true;
    }
}
void invalidateSnapshots(OptimizationPackage& package) {
    package.snapshotsValid = false;
}
double evaluatePathFromSnapshot(const vector<int>& path, const vector<SimSnapshot>& snapshots, int startIndex) {
    startIndex = max(0, min(startIndex, static_cast<int>(snapshots.size()) - 1));
    SimSnapshot snapshot = snapshots[startIndex];
    simulateUpgradePathFrom(path, snapshot, startIndex);
    return calculateScore(snapshot.resources);
}
double evaluatePath(const vector<int>& path, const SearchContext& context) {
    array<double, NUM_RESOURCES> testResources;
    array<int, NUM_UPGRADES + 1> testLevels;
    testResources = context.resources;
    testLevels = context.levels;
    simulateUpgradePath(path, testLevels, testResources);
    return calculateScore(testResources);
}
double calculateFinalPath(vector<int>& path, const array<int, NUM_UPGRADES + 1>& startingLevels = currentLevels){
    array<int, NUM_UPGRADES + 1>     simulationLevels(startingLevels);
    array<double, NUM_RESOURCES>  simulationResources(resourceCounts);

    simulateUpgradePath(path, simulationLevels, simulationResources, true);
    double simulationScore = calculateScore(simulationResources);
    printFormattedResults(path, simulationLevels, simulationResources, simulationScore);
    return simulationScore;
}
#if ENABLE_PARALLEL_SWAP
struct SwapSearchResult {
    bool found = false;
    long long order = 0;
    int indexA = -1;
    int indexB = -1;
    double score = 0;
};
long long pairsBeforeRow(int row, int searchLength) {
    return static_cast<long long>(row) * (2LL * searchLength - row - 1) / 2;
}
pair<int, int> pairFromOrder(long long order, int searchLength) {
    int low = 0;
    int high = searchLength - 2;

    while (low < high) {
        int mid = (low + high + 1) / 2;
        if (pairsBeforeRow(mid, searchLength) <= order) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    int first = low;
    int second = first + 1 + static_cast<int>(order - pairsBeforeRow(first, searchLength));
    return {first, second};
}
int chooseSwapThreadCount(long long candidateCount) {
    unsigned int hardwareThreads = thread::hardware_concurrency();
    if (candidateCount < PARALLEL_SWAP_MIN_CANDIDATES || hardwareThreads < 2) {
        return 1;
    }

    long long workSizedThreads = max(2LL, candidateCount / 5000);
    return static_cast<int>(min<long long>({hardwareThreads, PARALLEL_SWAP_MAX_THREADS, workSizedThreads}));
}
bool trySwapUpgradesParallel(OptimizationPackage& package, SearchContext& context, int startPos, int searchLength, long long candidateCount, int threadCount) {
    atomic<long long> earliestImprovement(candidateCount);
    vector<SwapSearchResult> results(threadCount);
    vector<thread> workers;
    workers.reserve(threadCount);

    for (int threadIndex = 0; threadIndex < threadCount; threadIndex++) {
        long long beginOrder = candidateCount * threadIndex / threadCount;
        long long endOrder = candidateCount * (threadIndex + 1) / threadCount;

        workers.emplace_back([&, threadIndex, beginOrder, endOrder]() {
            vector<int> localPath = package.path;
            auto [row, column] = pairFromOrder(beginOrder, searchLength);

            for (long long order = beginOrder; order < endOrder && order < earliestImprovement.load(memory_order_relaxed); order++) {
                int indexA = row + startPos;
                if (indexA >= searchLength) indexA -= searchLength;

                int indexB = column + startPos;
                if (indexB >= searchLength) indexB -= searchLength;

                if (localPath[indexA] != localPath[indexB]) {
                    swap(localPath[indexA], localPath[indexB]);
                    double testScore = evaluatePathFromSnapshot(localPath, package.snapshots, min(indexA, indexB));

                    if (testScore > package.score) {
                        results[threadIndex] = SwapSearchResult{true, order, indexA, indexB, testScore};

                        long long observed = earliestImprovement.load(memory_order_relaxed);
                        while (order < observed && !earliestImprovement.compare_exchange_weak(observed, order, memory_order_relaxed)) {}
                        break;
                    }

                    swap(localPath[indexA], localPath[indexB]);
                }

                column++;
                if (column >= searchLength) {
                    row++;
                    column = row + 1;
                }
            }
        });
    }

    for (thread& worker : workers) {
        worker.join();
    }

    SwapSearchResult bestResult;
    bestResult.order = candidateCount;
    for (const SwapSearchResult& result : results) {
        if (result.found && result.order < bestResult.order) {
            bestResult = result;
        }
    }

    if (!bestResult.found) {
        return false;
    }

    swap(package.path[bestResult.indexA], package.path[bestResult.indexB]);
    package.score = bestResult.score;
    context.logger.logImprovement("Swap", package.path, package.score);
    invalidateSnapshots(package);
    return true;
}
#endif
bool tryInsertUpgrade(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size();

    uniform_int_distribution<> positionDist(0, pathLength);
    int startPosition = positionDist(package.randomEngine);
    const int maxTypes = (NUM_RESOURCES * 2);
    ensureSnapshots(package, context);

    for (int i = 0; i < pathLength; i++) {
        int modulatedInsertPosition = i + startPosition;
        if (modulatedInsertPosition >= pathLength) modulatedInsertPosition -= pathLength;
        package.path.insert(package.path.begin() + modulatedInsertPosition, 0);

        int startingUpgradeType = package.randomEngine() % (NUM_RESOURCES * 2);
        for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {

            int modulatedUpgradeType = upgradeType + startingUpgradeType;
            if (modulatedUpgradeType >= maxTypes) modulatedUpgradeType -= maxTypes;
            package.path[modulatedInsertPosition] = modulatedUpgradeType;

            double testScore = evaluatePathFromSnapshot(package.path, package.snapshots, modulatedInsertPosition);

            if (testScore > package.score) {
                package.score = testScore;
                context.logger.logImprovement("Insert", package.path, package.score);
                invalidateSnapshots(package);
                return true;
            }
        }
        package.path.erase(package.path.begin() + modulatedInsertPosition);
    }
    package.deadInsert = true;
    return false;
}
bool tryRemoveUpgrade(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size() - 1;

    uniform_int_distribution<> swapDist(0, pathLength - 2);
    int startPos = swapDist(package.randomEngine);
    ensureSnapshots(package, context);

    for (int i = 0; i < pathLength; i++) {
        int removePos = i + startPos;
        if (removePos >= pathLength) removePos -= pathLength;
        int removedUpgrade = package.path[removePos];

        package.path.erase(package.path.begin() + removePos);
        double testScore = evaluatePathFromSnapshot(package.path, package.snapshots, removePos);

        if (testScore >= package.score) {
            package.score = testScore;
            context.logger.logImprovement("Remove", package.path, package.score);
            invalidateSnapshots(package);
            return true;
        }
        package.path.insert(package.path.begin() + removePos, removedUpgrade);
    }
    package.deadRemove = true;
    return false;

}
bool tryReplaceUpgrade(OptimizationPackage& package, SearchContext& context) {

    int pathLength = package.path.size();

    uniform_int_distribution<> positionDist(0, pathLength - 1); // -1 because we access index, not insert
    int startPosition = positionDist(package.randomEngine);
    const int maxTypes = (NUM_RESOURCES * 2);
    ensureSnapshots(package, context);

    for (int i = 0; i < pathLength; i++) {

        int targetIndex = i + startPosition;
        if (targetIndex >= pathLength) targetIndex -= pathLength;
        int originalType = package.path[targetIndex];

        int startingUpgradeType = package.randomEngine() % (NUM_RESOURCES * 2);
        
        for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {
            int modulatedUpgradeType = upgradeType + startingUpgradeType;
            if (modulatedUpgradeType >= maxTypes) modulatedUpgradeType -= maxTypes;

            if (modulatedUpgradeType == originalType) continue;

            package.path[targetIndex] = modulatedUpgradeType;

            double testScore = evaluatePathFromSnapshot(package.path, package.snapshots, targetIndex);

            if (testScore > package.score) {
                package.score = testScore;
                context.logger.logImprovement("Replace", package.path, package.score);
                invalidateSnapshots(package);
                return true;
            }
        }
        package.path[targetIndex] = originalType;
    }
    package.deadReplace = true;
    return false;
}
bool trySwapUpgrades(OptimizationPackage& package, SearchContext& context) {

    int pathLength = static_cast<int>(package.path.size()) - 1;
    int searchLength = pathLength - 1;
    if (searchLength < 2) {
        package.deadSwap = true;
        return false;
    }

    uniform_int_distribution<> swapDist(0, searchLength - 1);
    int startPos = swapDist(package.randomEngine);
    ensureSnapshots(package, context);

#if ENABLE_PARALLEL_SWAP
    long long candidateCount = static_cast<long long>(searchLength) * (searchLength - 1) / 2;
    int threadCount = chooseSwapThreadCount(candidateCount);
    if (threadCount > 1) {
        bool improved = trySwapUpgradesParallel(package, context, startPos, searchLength, candidateCount, threadCount);
        if (improved) return true;
        package.deadSwap = true;
        return false;
    }
#endif

    for (int i2 = 0; i2 < searchLength; i2++) { 
        for (int j2 = i2 + 1; j2 < searchLength; j2++) { 
            int i = i2 + startPos;
            if (i >= searchLength) i -= searchLength;
            int j = j2 + startPos;
            if (j >= searchLength) j -= searchLength;

            if (package.path[i] == package.path[j]) continue;
            swap(package.path[i], package.path[j]);

            double testScore = evaluatePathFromSnapshot(package.path, package.snapshots, min(i, j));

            if (testScore > package.score) {
                package.score = testScore;
                context.logger.logImprovement("Swap", package.path, package.score);
                invalidateSnapshots(package);
                return true;
            }
            swap(package.path[i], package.path[j]);
        }
    }
    package.deadSwap = true;
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
    ensureSnapshots(package, context);

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

        testScore = evaluatePathFromSnapshot(package.path, package.snapshots, rangeStart);
        if (testScore > package.score) {
            package.score = testScore;
            context.logger.logImprovement("Rotation", package.path, package.score);
            invalidateSnapshots(package);
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
    ensureSnapshots(package, context);

    for (int i = 0; i < maxIndex - 1; i++){
        int rangeStart = startOffset + i;
        if (rangeStart >= maxIndex - 1) rangeStart -= (maxIndex - 1);

        int remainingLength = maxIndex - (rangeStart + 2);
        uniform_int_distribution<> endOffsetDist(0, max(0, remainingLength)); 
        int endOffset = endOffsetDist(package.randomEngine);
        for (int j = 0; j < remainingLength; j++){ 
            int offsetJ = endOffset + j;
            if (offsetJ >= remainingLength + 1) offsetJ -= (remainingLength + 1);
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

                double testScore = evaluatePathFromSnapshot(package.path, package.snapshots, rangeStart);

                if (testScore > package.score) {
                    package.score = testScore;
                    context.logger.logImprovement("Rotation", package.path, package.score);
                    invalidateSnapshots(package);
                    return true;
                }
                //To undo Left(X), we Rotate Right(X)
                //To undo Right(X), we Rotate Left(X)
                if(rotateLeft)  rotate(rangeBeginIt, rangeEndIt - shiftAmount, rangeEndIt);
                else            rotate(rangeBeginIt, rangeBeginIt + shiftAmount, rangeEndIt);
            }
        }
    }
    package.deadRotate = true;
    return false;
}
void optimizeUpgradePath(OptimizationPackage& package, SearchContext& context, const int maxIterations = 10000) {

    random_device seed;
    mt19937 randomEngine(seed());

    int iterationCount = 0;
    int noImprovementStreak = 0;
    int InsertCount = 0, RemoveCount = 0, SwapCount = 0, RotationCount = 0, ReplaceCount = 0;
    double InsertRatio = 0, RemoveRatio = 0, SwapRatio = 0, RotationRatio = 0, ReplaceRatio = 0;
    auto benchmarkStart = chrono::high_resolution_clock::now();
    while (noImprovementStreak < maxIterations) {
        iterationCount++;
        bool improved = false;

        int strategy = iterationCount % 100;

        if(package.deadRotate){
            break;
        }
        else if (package.deadInsert && package.deadRemove && package.deadReplace && package.deadSwap)
                                                            {improved = tryRotateSubsequences(package, context); if (improved) RotationCount++;}
        else if (strategy < 15 && !package.deadInsert)      {improved = tryInsertUpgrade(package, context); if (improved) InsertCount++;}
        else if (strategy < 30 && !package.deadReplace)     {improved = tryReplaceUpgrade(package, context); if (improved) ReplaceCount++;}
        else if (strategy < 45 && !package.deadRemove)      {improved = tryRemoveUpgrade(package, context); if (improved) RemoveCount++;}
        else if (strategy < 47)                             {improved = tryRotateSubsequences(package, context); if (improved) RotationCount++;}
        else if (!package.deadSwap)                         {improved = trySwapUpgrades(package, context); if (improved) SwapCount++;}
        
        if (improved) {
            noImprovementStreak = 0;
            package.deadInsert = false;
            package.deadRemove = false;
            package.deadReplace = false;
            package.deadSwap = false;
            package.deadRotate = false;
            continue;
        } 
        else noImprovementStreak++;
    }
    //Data collection for analysis. Realistically needs its own manager.
    auto benchmarkEnd = chrono::high_resolution_clock::now();
    double elapsedSeconds = chrono::duration<double>(benchmarkEnd - benchmarkStart).count();
    int totalEvaluations = InsertCount + ReplaceCount + RemoveCount + SwapCount + RotationCount; 
    
    cout << "\n--- BENCHMARK RESULTS ---\n";
    cout << "Time elapsed: " << elapsedSeconds << " seconds\n";
    cout << "Total loop iterations: " << iterationCount << "\n";
    cout << "Throughput: " << (iterationCount / elapsedSeconds) << " iterations/second\n";
    cout << "-------------------------\n\n";
    // InsertRatio = (double)InsertCount / iterationCount;
    // ReplaceRatio = (double)ReplaceCount / iterationCount;
    // RemoveRatio = (double)RemoveCount / iterationCount;
    // SwapRatio = (double)SwapCount / iterationCount;
    // RotationRatio = (double)RotationCount / iterationCount;
    // cout << "Insert Count: " << InsertCount << "\n";
    // cout << "Replace Count: " << ReplaceCount << "\n";
    // cout << "Remove Count: " << RemoveCount << "\n";
    // cout << "Swap Count: " << SwapCount << "\n";
    // cout << "Rotation Count: " << RotationCount << "\n";
    // cout << "Insert ratio: " << InsertRatio << "\n";
    // cout << "Replace ratio: " << ReplaceRatio << "\n";
    // cout << "Remove ratio: " << RemoveRatio << "\n";
    // cout << "Swap ratio: " << SwapRatio << "\n";
    // cout << "Rotation ratio: " << RotationRatio << "\n";
}

int main() {
    upgradePath.reserve(500);
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

        array<int, NUM_UPGRADES + 1> activeLevels(currentLevels);
        if (isFullPath) {
            activeLevels = adjustFullPath(activeLevels);
        }

        double initial_score = calculateFinalPath(upgradePath, activeLevels);
        random_device seed;
        mt19937 randomEngine(seed());
        Logger logger(outputInterval);
        SearchContext context{logger, resourceCounts, activeLevels};
        OptimizationPackage package = {upgradePath, initial_score, move(randomEngine)};
        optimizeUpgradePath(package, context);
        upgradePath = move(package.path);
        
        finalScore = calculateFinalPath(upgradePath, activeLevels);
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
    array<int, NUM_UPGRADES + 1> activeLevels(currentLevels);
    if (isFullPath) {
        activeLevels = adjustFullPath(activeLevels);
    }

    double initial_score = calculateFinalPath(upgradePath, activeLevels);

    if (runOptimization) {
        random_device seed;
        mt19937 randomEngine(seed());
        Logger logger(outputInterval);
        SearchContext context{logger, resourceCounts, activeLevels};
        OptimizationPackage package = {upgradePath, initial_score, move(randomEngine)};
        optimizeUpgradePath(package, context);
        upgradePath = move(package.path);
    }
    finalScore = calculateFinalPath(upgradePath, activeLevels);
    
    return 0;
}   

