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
const int UNLOCKED_PETS = 52;
const int DLs = 1138;

vector<int> currentLevels = { 
    // Production Levels
    0, 1, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,   // Event Currency
    
    // Speed Levels
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,
    0, 0, 0,   // Event Currency
    
    0             // Dummy placeholder. Keep 0
};
vector<double> resourceCounts = { 
    0, 500000, 0,    
    0, 0, 0,  
    0/((500.0+DLs)/5.0), 0, 0,
    0, 0, 0/(UNLOCKED_PETS/100.0)                        
};

vector<int> upgradePath = {};

const bool isFullPath = false;
const bool runOptimization = true;   // Set false to only see the results and timings of your path
const bool endlessMode = false; //NOT WORKING YET-- Repeat optimization until *manually* stopped. Only save best result of all of them. Only works running locally

// END USER SETTINGS ---------------------------------------------------------------------
// ADVANCED SETTINGS ---------------------------------------------------------------------
const int outputInterval = 3000; //How often results may be printed (milliseconds)

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
constexpr int TOTAL_SECONDS = ((EVENT_DURATION_DAYS)*24*3600+(EVENT_DURATION_HOURS)*3600+(EVENT_DURATION_MINUTES)*60+EVENT_DURATION_SECONDS);

map<int, string> upgradeNames;
array<double, TOTAL_SECONDS> timeNeededSeconds{};

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
    static Proposal Replace(int index, int upgradeType, double score) {
        return Proposal{"Replace", score, index, 0, 0, upgradeType};
    }
};
void nameUpgrades() {
    for (int i = 0; i < NUM_RESOURCES; i++) {
        upgradeNames[i] = string(resourceNames[i]) + "_Level";
        upgradeNames[i + NUM_RESOURCES] = string(resourceNames[i]) + "_Speed";
    }
    upgradeNames[24] = "Complete";
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
double performUpgrade(vector<int>& levels, vector<double>& resources, int upgradeType, double& remainingTime) {

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
    
    int resourceType = upgradeType % NUM_RESOURCES;
    if (upgradeType == NUM_RESOURCES * 2) {
        resourceType = -1;
    }
    
    int newLevel = levels[upgradeType] + 1;
    double baseCost = newLevel * newLevel * newLevel * 100.0;

    if (upgradeType >= NUM_RESOURCES) baseCost *= 2.5;

    
    double cost[12] = {};
    double productionRates[12] = {};
    
    switch (resourceType) {
        case 0:
            cost[1] = baseCost * 2;
            break;
        case 1:
            cost[1] = baseCost;
            break;
        case 2:
            cost[1] = baseCost * 2;
            break;
        case 3:
            cost[0] = baseCost * 1.5;
            cost[1] = baseCost * 1.5;
            break;
        case 4:
            cost[1] = baseCost * 2;
            break;
        case 5:
            cost[1] = baseCost * 1.5;
            cost[2] = baseCost * 1.5;
            break;
        case 6:
            cost[3] = baseCost * 5;
            break;
        case 7:
            cost[1] = baseCost * 1.5;
            cost[4] = baseCost * 1.5;
            break;
        case 8:
            cost[5] = baseCost * 5;
            break;
        case 9:
            cost[3] = baseCost * 2.5;
            cost[4] = baseCost * 2.5;
            break;
        case 10:
            cost[3] = baseCost * 3;
            cost[4] = baseCost * 3;
            cost[5] = baseCost * 3;
            break;
        case 11:
            cost[5] = baseCost * 2.5;
            cost[7] = baseCost * 2.5;
            break;
    }
    
    productionRates[0] = levels[0] * cycleTimeMultiplier[0] * speedMultipliers[levels[12]]; // Turning this into a loop worsenes performance
    productionRates[1] = levels[1] * cycleTimeMultiplier[1] * speedMultipliers[levels[13]]; 
    productionRates[2] = levels[2] * cycleTimeMultiplier[2] * speedMultipliers[levels[14]];
    productionRates[3] = levels[3] * cycleTimeMultiplier[3] * speedMultipliers[levels[15]];
    productionRates[4] = levels[4] * cycleTimeMultiplier[4] * speedMultipliers[levels[16]];
    productionRates[5] = levels[5] * cycleTimeMultiplier[5] * speedMultipliers[levels[17]];
    productionRates[6] = levels[6] * cycleTimeMultiplier[6] * speedMultipliers[levels[18]];
    productionRates[7] = levels[7] * cycleTimeMultiplier[7] * speedMultipliers[levels[19]];
    productionRates[8] = levels[8] * cycleTimeMultiplier[8] * speedMultipliers[levels[20]];
    productionRates[9] = levels[9] * cycleTimeMultiplier[9] * speedMultipliers[levels[21]];
    productionRates[10] = levels[10] * cycleTimeMultiplier[10] * speedMultipliers[levels[22]];
    productionRates[11] = levels[11] * cycleTimeMultiplier[11] * speedMultipliers[levels[23]];
    
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
        resources[10] += productionRates[10] * timeNeeded;
        resources[11] += productionRates[11] * timeNeeded;
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
    resources[10] += productionRates[10] * timeNeeded - cost[10];
    resources[11] += productionRates[11] * timeNeeded - cost[11];

    levels[upgradeType]++;
    
    return timeNeeded;
}
double simulateUpgradePath(vector<int>& path, vector<int>& levels, vector<double>& resources, bool display = false) {
    thread_local double time;
    time = TOTAL_SECONDS;
    for (auto upgradeType : path) {
        if (time < 1e-3) return 0;
        if (upgradeType >= NUM_RESOURCES && levels[upgradeType] >= 10) return INFINITY_VALUE;
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
bool tryInsertUpgrade(OptimizationPackage& package, SearchContext& context, Proposal* outProposal = nullptr) {

    int pathLength = package.path.size();
    thread_local vector<int> candidatePath;

    uniform_int_distribution<> positionDist(0, pathLength);
    int startPosition = positionDist(package.randomEngine);
    const int maxTypes = (NUM_RESOURCES * 2);

    for (int i = 0; i < pathLength; i++) {
        candidatePath = package.path;
        int modulatedInsertPosition = (i + startPosition) % pathLength;
        candidatePath.insert(candidatePath.begin() + modulatedInsertPosition, 0);

        int startingUpgradeType = package.randomEngine() % (NUM_RESOURCES * 2);
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
bool tryReplaceUpgrade(OptimizationPackage& package, SearchContext& context, Proposal* outProposal = nullptr) {

    int pathLength = package.path.size();
    thread_local vector<int> candidatePath;

    uniform_int_distribution<> positionDist(0, pathLength - 1); // -1 because we access index, not insert
    int startPosition = positionDist(package.randomEngine);
    const int maxTypes = (NUM_RESOURCES * 2);

    for (int i = 0; i < pathLength; i++) {
        // We do not resize the vector here
        candidatePath = package.path; 
        int targetIndex = (i + startPosition) % pathLength;
        
        // Store original to ensure we are actually changing it
        int originalType = candidatePath[targetIndex];

        int startingUpgradeType = package.randomEngine() % (NUM_RESOURCES * 2);
        
        for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {
            int modulatedUpgradeType = (upgradeType + startingUpgradeType) % maxTypes;

            // Don't replace an upgrade with itself
            if (modulatedUpgradeType == originalType) continue;

            // PERFORM REPLACEMENT
            candidatePath[targetIndex] = modulatedUpgradeType;

            double testScore = evaluatePath(candidatePath, context);

            if (testScore > package.score) {
                if (outProposal) *outProposal = Proposal::Replace(targetIndex, modulatedUpgradeType, testScore);
                
                // Commit change to package
                package.path = candidatePath; 
                package.score = testScore;
                context.logger.logImprovement("Replace", package.path, package.score);
                return true;
            }
        }
    }
    package.deadMoves.insert("Replace");
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
    nameUpgrades();
    preprocessBusyTimes(busyTimesStart, busyTimesEnd);

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
    finalScore =calculateFinalPath(upgradePath);
    
    return 0;
}   

