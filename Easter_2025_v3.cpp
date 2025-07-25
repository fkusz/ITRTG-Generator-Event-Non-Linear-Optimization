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
// (WIP) added the option for absolute timestamps for upgrade buy-times rather than timestamps relative to the start of the sim
// (WIP) add parallelization for a isngle solution to produce a near-instant single solution
// (WIP) rare bug where the score can become unbounded and skyrocket

#include <windows.h>
#include <iostream>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <random> 
#include <array>
#include <stop_token>
#include <thread>
#include <atomic>
#include <mutex>
#include <unordered_set>
#include <shared_mutex>
using namespace std;
typedef long long ll;

// USER SETTINGS ---------------------------------------------------------------------
const int EVENT_DURATION_DAYS = 14;
const int EVENT_DURATION_HOURS = 0;
const int EVENT_DURATION_MINUTES = 0;
const int EVENT_DURATION_SECONDS = 0;
const int UNLOCKED_PETS = 37;
const int DLs = 369;

vector<int> currentLevels = { 
    // Production Levels
    0, 1, 0,
    0, 0, 0,
    0, 0, 0,
    0,            // Event Currency
    
    // Speed Levels
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0,            // Event Currency
    
    0             // Dummy placeholder
};
vector<double> resourceCounts = { 
    0,  500000, 0,    
    0,  0, 0,   
    0,   
    0/((500.0+DLs)/5.0), 
    0/(UNLOCKED_PETS/100.0),                        
    0,                        
};

vector<int> upgradePath = {};

const bool isFullPath = false;
const bool allowSpeedUpgrades = true; 
const bool runOptimization = true;   

// END USER SETTINGS ---------------------------------------------------------------------
// ADVANCED SETTINGS ---------------------------------------------------------------------
const int outputInterval = 1000; //How often results may be printed (milliseconds)
const bool multithreading = false; // If true, run on multiple CPU cores. (Probably not worth it yet)

const double EVENT_CURRENCY_WEIGHT = 100e-5; 
const double FREE_EXP_WEIGHT = 6e-5;     
const double PET_STONES_WEIGHT = 3.5e-5;  
const double GROWTH_WEIGHT = 8e-5;       

const vector<double> busyTimesStart =   {}; // The first hours, relative to the start of the simulation, you'll be unable to purchase upgrades
const vector<double> busyTimesEnd =     {}; // The first hours, relative to the start of the simulation, you'll be able to purchase upgrades again
// END ADVANCED SETTINGS ------------------------------------------------------------------
// PROGRAM SETTINGS -----------------------------------------------------------------------
constexpr array<const char*, 10> resourceNames = {"Chair", "Bucket", "Goggles", "Water_Gun", "Surfboard", "Sunglasses", "PET_STONES", "FREE_EXP", "GROWTH", "EVENT_CURRENCY"};
constexpr double INFINITY_VALUE = (1e100);
constexpr int NUM_RESOURCES = resourceNames.size();
constexpr int TOTAL_SECONDS = ((EVENT_DURATION_DAYS)*24*3600+(EVENT_DURATION_HOURS)*3600+(EVENT_DURATION_MINUTES)*60+EVENT_DURATION_SECONDS);
int num_cores = thread::hardware_concurrency()/2;

map<int, string> upgradeNames;
array<double, TOTAL_SECONDS> timeNeededSeconds{};

// END PROGRAM SETTINGS -----------------------------------------------------------------
// UTILITY FUNCTIONS --------------------------------------------------------------------
template <typename T>
void printVector(vector<T>& x, ostream& out = cout) {
    for (auto item : x) out << item << ",";
    out << "\n";
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
            out << "Improved path: (" << type << ") ";
            printVector(path, out);
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
struct Proposal {
    string type;
    double newScore;
    int indexA = 0;
    int indexB = 0;
    int rotateIndex = 0;
    int upgrade = 0;

    static Proposal Insert(int index, int upgradeType, double score) {
        return Proposal{"Insert", score, index, 0, 0, upgradeType};
    }
    static Proposal Remove(int index,  double score) {
        return Proposal{"Remove", score, index, 0, 0, 0};
    }
    static Proposal Swap(int indexA, int IndexB, double score) {
        return Proposal{"Swap", score, indexA, IndexB, 0, 0};
    }
    static Proposal Rotate(int indexA, int indexB, int rotateIndex, double score) {
        return Proposal{"Rotate", score, indexA, indexB, rotateIndex, 0};
    }

    // Add other helpers as needed
};
void nameUpgrades() {
    for (int i = 0; i < NUM_RESOURCES; i++) {
        upgradeNames[i] = string(resourceNames[i]) + "_Level";
        upgradeNames[i + NUM_RESOURCES] = string(resourceNames[i]) + "_Speed";
    }
    upgradeNames[20] = "Complete";
}
vector <int> adjustFullPath(vector<int>& levels){ // This is not thread-safe
    levels[1]--; // Adjust first level because it always starts at 1
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
    cout << "Upgrade Path: ";
    printVector(path);
    cout << "Final Resource Counts: ";
    printVector(simulationResources);
    cout << "Final Upgrade Levels: ";
    printVector(simulationLevels);

    cout << "Event Currency: " << simulationResources[9] << "\n";
    cout << "Free Exp (" << DLs << " DLs): " 
        << simulationResources[7] * (500.0 + DLs) / 5.0 
        << " (" << simulationResources[7] << " levels * cycles)" << "\n";
    cout << "Pet Stones: " << simulationResources[6] << "\n";
    cout << "Growth (" << UNLOCKED_PETS << " pets): " 
        << simulationResources[8] * UNLOCKED_PETS / 100.0 
        << " (" << simulationResources[8] << " levels * cycles)" << "\n";
    cout << "Score: " << finalScore << "\n";
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
        randomPath.push_back(rand() % NUM_RESOURCES + 10 * (rand() % 2));
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
void setThreadName(const std::string& name) {
    HRESULT hr = SetThreadDescription(GetCurrentThread(), std::wstring(name.begin(), name.end()).c_str());
    if (FAILED(hr)) {
        std::cerr << "Failed to set thread name: " << std::hex << hr << "\n";
    }
}
// END UTILITY FUNCTIONS ----------------------------------------------------------------
// ALGORITHM FUNCTIONS ------------------------------------------------------------------
double performUpgrade(vector<int>& levels, vector<double>& resources, int upgradeType, double& remainingTime) {

    constexpr double cycleTimeMultiplier[10] = {
        1.0/3.0,    
        1.0,        
        1.0/3.0,    
        1.0/3.0,    
        1.0/3.0,    
        1.0/3.0,    
        1.0/1200.0, 
        1.0/2500.0, 
        1.0/1800.0, 
        1.0/5000.0  
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
    
    int resourceType = upgradeType % NUM_RESOURCES;
    if (upgradeType == NUM_RESOURCES * 2) {
        resourceType = -1;
    }
    
    int newLevel = levels[upgradeType] + 1;
    double baseCost = (3.0 * newLevel * newLevel * newLevel + 1.0) * 100.0;

    if (upgradeType >= NUM_RESOURCES) baseCost *= 2.0;

    
    double cost[10] = {};
    double productionRates[10] = {};
    
    switch (resourceType) {
        case 0:
            cost[1] = baseCost * 10;
            break;
        case 1:
            cost[1] = baseCost * 0.8;
            break;
        case 2:
            cost[1] = baseCost;
            break;
        case 3:
            cost[0] = baseCost;
            break;
        case 4:
            cost[0] = baseCost;
            break;
        case 5:
            cost[1] = baseCost;
            cost[2] = baseCost;
            break;
        case 6:
            cost[0] = baseCost * 0.7;
            cost[3] = baseCost * 0.5;      
            break;
        case 7:
            cost[0] = baseCost;    
            cost[4] = baseCost * 3;
            break;
        case 8:
            cost[2] = baseCost * 1.2;
            cost[5] = baseCost;
            break;
        case 9:
            cost[3] = baseCost;
            cost[4] = baseCost;
            cost[5] = baseCost;
            break;
    }
    
    productionRates[0] = levels[0] * cycleTimeMultiplier[0] * speedMultipliers[levels[10]]; // Turning this into a loop worsenes performance
    productionRates[1] = levels[1] * cycleTimeMultiplier[1] * speedMultipliers[levels[11]]; 
    productionRates[2] = levels[2] * cycleTimeMultiplier[2] * speedMultipliers[levels[12]];
    productionRates[3] = levels[3] * cycleTimeMultiplier[3] * speedMultipliers[levels[13]];
    productionRates[4] = levels[4] * cycleTimeMultiplier[4] * speedMultipliers[levels[14]];
    productionRates[5] = levels[5] * cycleTimeMultiplier[5] * speedMultipliers[levels[15]];
    productionRates[6] = levels[6] * cycleTimeMultiplier[6] * speedMultipliers[levels[16]];
    productionRates[7] = levels[7] * cycleTimeMultiplier[7] * speedMultipliers[levels[17]];
    productionRates[8] = levels[8] * cycleTimeMultiplier[8] * speedMultipliers[levels[18]];
    productionRates[9] = levels[9] * cycleTimeMultiplier[9] * speedMultipliers[levels[19]];
    
    // Calculate time needed for the upgrade
    double timeNeeded = 0;
    double timeElapsed = TOTAL_SECONDS - remainingTime;
    for (int i = 0; i < NUM_RESOURCES; i++) {
        double neededResources = cost[i] - resources[i];
        if (neededResources <= 0) continue;
        if (productionRates[i] == 0) {
            timeNeeded = INFINITY_VALUE; 
            break;
        }
        timeNeeded = max(timeNeeded, neededResources / productionRates[i]);
    }

    int busyLookupIndex = static_cast<int>(timeElapsed + timeNeeded);
    if (0 <= busyLookupIndex && busyLookupIndex < TOTAL_SECONDS){
        timeNeeded += timeNeededSeconds[busyLookupIndex];
    };

    if (timeNeeded >= remainingTime || upgradeType == (2 * NUM_RESOURCES)) {
        timeNeeded = remainingTime;
        resources[0] += productionRates[0] * timeNeeded;
        resources[1] += productionRates[1] * timeNeeded; 
        resources[2] += productionRates[2] * timeNeeded;
        resources[3] += productionRates[3] * timeNeeded;
        resources[4] += productionRates[4] * timeNeeded;
        resources[5] += productionRates[5] * timeNeeded;
        resources[6] += productionRates[6] * timeNeeded;
        resources[7] += productionRates[7] * timeNeeded;
        resources[8] += productionRates[8] * timeNeeded;
        resources[9] += productionRates[9] * timeNeeded;
        return timeNeeded;
    }

    resources[0] += productionRates[0] * timeNeeded - cost[0];
    resources[1] += productionRates[1] * timeNeeded - cost[1];
    resources[2] += productionRates[2] * timeNeeded - cost[2];
    resources[3] += productionRates[3] * timeNeeded - cost[3];
    resources[4] += productionRates[4] * timeNeeded - cost[4];
    resources[5] += productionRates[5] * timeNeeded - cost[5];
    resources[6] += productionRates[6] * timeNeeded - cost[6];
    resources[7] += productionRates[7] * timeNeeded - cost[7];
    resources[8] += productionRates[8] * timeNeeded - cost[8];
    resources[9] += productionRates[9] * timeNeeded - cost[9];

    levels[upgradeType]++;
    
    return timeNeeded;
}
double simulateUpgradePath(vector<int>& path, vector<int>& levels, vector<double>& resources, bool display = false) {
    thread_local double time;
    time = TOTAL_SECONDS;
    for (auto upgradeType : path) {
        if (time < 1e-3) return 0;
        
        double timeTaken = performUpgrade(levels, resources, upgradeType, time);
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

    score += (min(resources[9], 10000.0) + max(0.0, (resources[9] - 10000)) * 0.01) * (EVENT_CURRENCY_WEIGHT);
    
    score += resources[7] * (FREE_EXP_WEIGHT);     // Free EXP
    score += resources[8] * (GROWTH_WEIGHT);   // Growth
    score += resources[6] * (PET_STONES_WEIGHT);       // Pet Stones
    
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
void calculateFinalPath(vector<int>& path){
    vector<int>     simulationLevels(currentLevels);
    vector<double>  simulationResources(resourceCounts);

    simulateUpgradePath(path, simulationLevels, simulationResources, true);
    double simulationScore = calculateScore(simulationResources);
    printFormattedResults(path, simulationLevels, simulationResources, simulationScore);
}
bool tryInsertUpgrade(OptimizationPackage& package, SearchContext& context, Proposal* outProposal = nullptr) {

    int pathLength = package.path.size();
    thread_local vector<int> candidatePath;

    uniform_int_distribution<> positionDist(0, pathLength);
    int startPosition = positionDist(package.randomEngine);
    const int maxTypes = (allowSpeedUpgrades ? NUM_RESOURCES * 2 : NUM_RESOURCES);

    for (int i = 0; i < pathLength; i++) {
        candidatePath = package.path;
        int modulatedInsertPosition = (i + startPosition) % pathLength;
        candidatePath.insert(candidatePath.begin() + modulatedInsertPosition, 0);

        int startingUpgradeType = package.randomEngine() % (allowSpeedUpgrades ? NUM_RESOURCES * 2 : NUM_RESOURCES);
        for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {

            int modulatedUpgradeType = (upgradeType + startingUpgradeType) % maxTypes;
            candidatePath[modulatedInsertPosition] = modulatedUpgradeType;

            double testScore = evaluatePath(candidatePath, context);

            if (testScore > package.score) {
                if (outProposal) *outProposal = Proposal::Insert(modulatedInsertPosition, modulatedUpgradeType, testScore);
                package.path.insert(package.path.begin() + modulatedInsertPosition, modulatedUpgradeType);
                package.score = testScore;
                context.logger.logImprovement("Insert", package.path, package.score);
                return true;
            }
        }
    }
    package.deadMoves.insert("Insert");
    return false;
}
bool tryRemoveUpgrade(OptimizationPackage& package, SearchContext& context, Proposal* outProposal = nullptr) {

    int pathLength = package.path.size() - 1;
    thread_local vector<int> candidatePath;

    uniform_int_distribution<> swapDist(0, pathLength - 2);
    int startPos = swapDist(package.randomEngine);

    for (int i = 0; i < pathLength; i++) {
        int removePos = (i + startPos) % (pathLength);

        if (!allowSpeedUpgrades && package.path[removePos] >= NUM_RESOURCES) continue;

        candidatePath = package.path;
        candidatePath.erase(candidatePath.begin() + removePos);
        double testScore = evaluatePath(candidatePath, context);

        if (testScore >= package.score) {
            if (outProposal) *outProposal = Proposal::Remove(removePos, testScore);
            package.score = testScore;
            package.path = candidatePath;
            context.logger.logImprovement("Remove", package.path, package.score);
            return true;
        }
    }
    package.deadMoves.insert("Remove");
    return false;

}
bool trySwapUpgrades(OptimizationPackage& package, SearchContext& context, Proposal* outProposal = nullptr) {

    int pathLength = package.path.size() - 1;
    thread_local vector<int> candidatePath;
    candidatePath = package.path;
    double testScore;

    uniform_int_distribution<> swapDist(0, pathLength - 2);
    int startPos = swapDist(package.randomEngine);

    for (int i2 = 0; i2 < pathLength - 1; i2++) { 
        for (int j2 = i2 + 1; j2 < pathLength - 1; j2++) { 
            int i = (i2 + startPos) % (pathLength - 1);
            int j = (j2 + startPos) % (pathLength - 1);

            if (candidatePath[i] == candidatePath[j]) continue;
            swap(candidatePath[i], candidatePath[j]);

            testScore = evaluatePath(candidatePath, context);

            if (testScore > package.score) {
                if (outProposal) *outProposal = Proposal::Swap(i, j, testScore);
                package.path = candidatePath;
                package.score = testScore;
                context.logger.logImprovement("Swap", package.path, package.score);
                return true;
            }

            swap(candidatePath[i], candidatePath[j]);
        }
    }
    package.deadMoves.insert("Swap");
    return false;
}
bool tryRotateSubsequences(OptimizationPackage& package, SearchContext& context, Proposal* outProposal = nullptr) {

    int pathLength = package.path.size() - 1;
    thread_local vector<int> candidatePath;
    double testScore;

    uniform_int_distribution<> rotateDist(0, pathLength - 3);
    int i = rotateDist(package.randomEngine);
    uniform_int_distribution<> rotateDist2(i+2, pathLength - 1);
    int j = rotateDist2(package.randomEngine);

    for (int k = 0; k < j-i; k++) {
        candidatePath = package.path;

        int offset = (k + 2) / 2;
        bool isLeft = (k % 2 == 0); 

        if(isLeft)  rotate(candidatePath.begin() + i, candidatePath.begin() + i + offset, candidatePath.begin() + j + 1); // Left rotation
        else        rotate(candidatePath.begin() + i, candidatePath.begin() + j - offset + 1, candidatePath.begin() + j + 1); // Right rotation

        testScore = evaluatePath(candidatePath, context);
        int rotationPos = isLeft ? i + offset: j - offset + 1;
        if (testScore > package.score) {
            if (outProposal) *outProposal = Proposal::Rotate(i, j + 1, rotationPos, testScore);
            package.path = candidatePath;
            package.score = testScore;
            context.logger.logImprovement("Rotation", package.path, package.score);
            return true;
        }
    }
    return false;
}
bool exhaustRotateSubsequences(OptimizationPackage& package, SearchContext& context, Proposal* outProposal = nullptr) {


    int pathLength = package.path.size() - 1;
    int maxIndex = pathLength - 1;
    thread_local vector<int> candidatePath;
    double testScore;

    uniform_int_distribution<> rotateDist(0, maxIndex - 2);
    int i = rotateDist(package.randomEngine);

    for (int i2 = 0; i2 < maxIndex - 1; i2++){
        int i3 = (i + i2) % (maxIndex - 1);
        uniform_int_distribution<> rotateDist2(0, maxIndex-i3); 
        int j = rotateDist2(package.randomEngine);
        for (int j2 = 0; j2 < maxIndex - i3 - 1; j2++){ 
            int j3 = i3 + 2 + ((j + j2) % (maxIndex - i3 - 1));
            for (int k = 0; k < j3-i3; k++) {
                candidatePath = package.path;

                int offset = (k + 2) / 2;
                bool isLeft = (k % 2 == 0);

                if(isLeft)  rotate(candidatePath.begin() + i3, candidatePath.begin() + i3 + offset, candidatePath.begin() + j3 + 1);
                else        rotate(candidatePath.begin() + i3, candidatePath.begin() + j3 - offset + 1, candidatePath.begin() + j3 + 1);

                double testScore = evaluatePath(candidatePath, context);
                int rotationPos = isLeft ? i3 + offset: j3 - offset + 1;
                if (testScore > package.score) {
                    if (outProposal) *outProposal = Proposal::Rotate(i3, j3+1, rotationPos, testScore);
                    package.path = candidatePath;
                    package.score = testScore;
                    context.logger.logImprovement("Rotation", package.path, package.score);
                    return true;
                }
            }

        }
    }
    package.deadMoves.insert("Rotate");
    return false;
}
void optimizeUpgradePath(OptimizationPackage& package, SearchContext& context, const int maxIterations = 1000) {

    random_device seed;
    mt19937 randomEngine(seed());

    int iterationCount = 0;
    int noImprovementStreak = 0;
    int InsertCount = 0, RemoveCount = 0, SwapCount = 0, RotationCount = 0;
    double InsertRatio = 0, RemoveRatio = 0, SwapRatio = 0, RotationRatio = 0;
    while (noImprovementStreak < maxIterations) {
        iterationCount++;
        bool improved = false;

        int strategy = iterationCount % 100;

        if(package.deadMoves.contains("Rotation")){
            break;
        }
        else if (package.deadMoves.contains("Insert") && package.deadMoves.contains("Remove") && package.deadMoves.contains("Swap"))
                                                                            {improved = tryRotateSubsequences(package, context); if (improved) RotationCount++;}
        else if (strategy < 15 && !package.deadMoves.contains("Insert"))    {improved = tryInsertUpgrade(package, context); if (improved) InsertCount++;}
        else if (strategy < 30 && !package.deadMoves.contains("Remove"))    {improved = tryRemoveUpgrade(package, context); if (improved) RemoveCount++;}
        else if (strategy < 32)                                             {improved = tryRotateSubsequences(package, context); if (improved) RotationCount++;}
        else if (!package.deadMoves.contains("Swap"))                       {improved = trySwapUpgrades(package, context); if (improved) SwapCount++;}
        
        if (improved) {
            noImprovementStreak = 0;
            package.deadMoves.clear();
            continue;
        } 
        else noImprovementStreak++;
    }
    // Data collection for analysis. Realistically needs its own manager.
    // InsertRatio = (double)InsertCount / iterationCount;
    // RemoveRatio = (double)RemoveCount / iterationCount;
    // SwapRatio = (double)SwapCount / iterationCount;
    // RotationRatio = (double)RotationCount / iterationCount;
    // cout << "Insert Count: " << InsertCount << "\n";
    // cout << "Remove Count: " << RemoveCount << "\n";
    // cout << "Swap Count: " << SwapCount << "\n";
    // cout << "Rotation Count: " << RotationCount << "\n";
    // cout << "Insert ratio: " << InsertRatio << "\n";
    // cout << "Remove ratio: " << RemoveRatio << "\n";
    // cout << "Swap ratio: " << SwapRatio << "\n";
    // cout << "Rotation ratio: " << RotationRatio << "\n";
}
vector<int> parallelOptimization() {
    random_device rd;
    vector<thread> threads;
    vector<OptimizationPackage> results(num_cores);

    cout << "Starting optimization on " << num_cores << " independent paths..." << "\n";

    for (int i = 0; i < num_cores; ++i) {
        threads.emplace_back([&, i]() {
            vector<int> threadPath;
            threadPath.reserve(500);
            threadPath = generateRandomPath();
            mt19937 threadRandomEngine(rd() ^ (hash<thread::id>{}(this_thread::get_id()) + i));
            Logger threadLogger(outputInterval);
            SearchContext threadContext{threadLogger, resourceCounts, currentLevels};
            OptimizationPackage threadPackage = {move(threadPath), 0, move(threadRandomEngine)};
            optimizeUpgradePath(threadPackage, threadContext);
            results[i] = move(threadPackage);
        });
    }

    for (auto& t : threads) {
        t.join();
    }
    for (auto& t : results) {
        cout << "Thread score: " << t.score << "\n";
    }

    // Pick the best result
    auto bestResult = max_element(results.begin(), results.end(), [](const auto& a, const auto& b) {
        return a.score < b.score;
    });

    return bestResult->path;
}
vector<int> masterWorkerOptimization() {
    mutex proposalMutex;
    mutex deadMoveMutex;
    shared_mutex packageMutex;
    atomic<uint64_t> currentVersion = 0; 
    optional<Proposal> sharedProposal;
    optional<string> sharedDeadMove;
    vector<int> masterPath;
    masterPath.reserve(500);
    masterPath = generateRandomPath();
    OptimizationPackage masterPackage = {masterPath, 0, mt19937(random_device{}())};
    setThreadName("Master");
    Logger masterLogger(outputInterval);

    // Worker thread logic
    auto workerFunc = [&](stop_token st, const string defaultStrategy) {
        mt19937 rng(random_device{}() ^ hash<thread::id>{}(this_thread::get_id()));
        setThreadName("WorkerThread");
        NullStream nullStream;
        Logger noopLogger(0, nullStream);
        SearchContext context(noopLogger, resourceCounts, currentLevels);
        Proposal proposal;
        string strategy;
        auto threadID = this_thread::get_id();
        OptimizationPackage localPackage;
        while (!st.stop_requested()) {
            bool improved = false;
            strategy = defaultStrategy;
            uint64_t myVersion = currentVersion.load();
            {
            shared_lock packageLock(packageMutex);
            localPackage = masterPackage;
            }

            if (!localPackage.deadMoves.contains(strategy));
            else if (!localPackage.deadMoves.contains("Insert")) strategy = "Insert";
            else if (!localPackage.deadMoves.contains("Swap")) strategy = "Swap";
            else if (!localPackage.deadMoves.contains("Remove")) strategy = "Remove";
            else if (!localPackage.deadMoves.contains("Rotate")) strategy = "Rotate";

            if      (strategy == "Insert") improved = tryInsertUpgrade(localPackage, context, &proposal);
            else if (strategy == "Swap")   improved = trySwapUpgrades(localPackage, context, &proposal);
            else if (strategy == "Remove") improved = tryRemoveUpgrade(localPackage, context, &proposal);
            else if (strategy == "Rotate") improved = exhaustRotateSubsequences(localPackage, context, &proposal);

            if (improved) {
                lock_guard proposalLock(proposalMutex);
                if (!sharedProposal.has_value()) {
                    if (myVersion != currentVersion.load()) continue;
                    sharedProposal = proposal;
                    //cout << "Worker Thread " << threadID << " found proposal " << strategy << "\n";
                    this_thread::sleep_for(1ms);
                }
            }
            else{
                lock_guard deadMoveLock(deadMoveMutex);
                if (!sharedDeadMove.has_value()) {
                    if (myVersion != currentVersion.load()) continue;
                    sharedDeadMove = strategy;
                    //cout << "Worker Thread " << threadID << " found dead move " << strategy << "\n";
                    this_thread::sleep_for(1ms);
                }
            }
        }
    };

    // Launch worker threads
    vector<jthread> workers;
    vector<string> strategies = {"Insert", "Remove", "Swap", "Swap", "Swap", "Swap", "Swap", "Swap", "Swap", "Swap", "Swap", "Swap", "Swap", "Swap", "Rotate"};
    for (int i = 0; i < num_cores - 1; ++i) {
        string preferredStrategy = strategies[i];
        workers.emplace_back([&, preferredStrategy](stop_token st) {workerFunc(st, preferredStrategy);});
    }

    while (true) {
        if (masterPackage.deadMoves.contains("Rotate")) break;
        lock_guard proposalLock(proposalMutex);
        if (sharedProposal.has_value()) {
            lock_guard deadMoveLock(deadMoveMutex);
            unique_lock packageLock(packageMutex);
            if (sharedProposal.value().type == "Insert") masterPackage.path.insert(masterPackage.path.begin() + sharedProposal.value().indexA, sharedProposal.value().upgrade);
            else if (sharedProposal.value().type == "Swap") swap(masterPackage.path[sharedProposal.value().indexA],masterPackage.path[sharedProposal.value().indexB]);
            else if (sharedProposal.value().type == "Remove") masterPackage.path.erase(masterPackage.path.begin() + sharedProposal.value().indexA);
            else if (sharedProposal.value().type == "Rotate") rotate(masterPackage.path.begin() + sharedProposal.value().indexA, masterPackage.path.begin() + sharedProposal.value().rotateIndex, masterPackage.path.begin() + sharedProposal.value().indexB);
            
            masterPackage.score = sharedProposal.value().newScore;
            masterPackage.deadMoves.clear();
            masterLogger.logImprovement(sharedProposal.value().type, masterPackage.path, masterPackage.score);
            sharedProposal.reset();
            sharedDeadMove.reset();
            currentVersion++;
            continue;
        }
        lock_guard deadMoveLock(deadMoveMutex);
        if (sharedDeadMove.has_value()) {
            unique_lock packageLock(packageMutex);
            masterPackage.deadMoves.insert(sharedDeadMove.value());
            sharedDeadMove.reset();
            currentVersion++;
            continue;
        }
    };

    return masterPackage.path;
}
int main() {
    upgradePath.reserve(500);
    vector<int> levelsCopy(currentLevels);
    
    nameUpgrades();
    preprocessBusyTimes(busyTimesStart, busyTimesEnd);

    if (runOptimization) {
        if (multithreading && num_cores > 2) {
            upgradePath = parallelOptimization();
        }
        else {

            // upgradePath = masterWorkerOptimization(); // Unfinished Speed Optimizations
            random_device seed;
            mt19937 randomEngine(seed());
            Logger logger(outputInterval);
            SearchContext context{logger, resourceCounts, currentLevels};
            OptimizationPackage package = {generateRandomPath(), 0, move(randomEngine)};
            optimizeUpgradePath(package, context);
            upgradePath = move(package.path);
        }
    }
    
    if (upgradePath.empty()) {          
        upgradePath = generateRandomPath();
    }
    if (isFullPath) {
        levelsCopy = adjustFullPath(levelsCopy);
    }

    calculateFinalPath(upgradePath);
    return 0;
}   

