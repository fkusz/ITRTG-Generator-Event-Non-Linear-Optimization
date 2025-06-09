#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <random>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <ctime> 
#include <array>
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
    0/(UNLOCKED_PETS/100.0),   
    0/((500.0+DLs)/5.0), 
    0,                        
    0,                        
};

vector<int> upgradePath = {};

const bool isFullPath = false;
const bool allowSpeedUpgrades = true; 
const bool runOptimization = true;   

const int outputInterval = 10; 
// END USER SETTINGS ---------------------------------------------------------------------
// ADVANCED SETTINGS ---------------------------------------------------------------------
const double EVENT_CURRENCY_WEIGHT = 100e-5; 
const double FREE_EXP_WEIGHT = 6e-5;     
const double PET_STONES_WEIGHT = 3.5e-5;  
const double GROWTH_WEIGHT = 8e-5;       

const vector<double> busyTimesStart = {
    16,40,64,88,112,136,160,184,208,232,256,280,304,328
}; 
const vector<double> busyTimesEnd = {
    24,48,72,96,120,144,168,192,216,240,264,288,312,336
};
// END ADVANCED SETTINGS ------------------------------------------------------------------
// PROGRAM SETTINGS -----------------------------------------------------------------------
constexpr std::array<const char*, 10> resourceNames = {"Red", "White", "Blue", "Green", "Colorful", "Yellow", "GROWTH", "FREE_EXP", "PET_STONES", "EVENT_CURRENCY"};
constexpr double INFINITY_VALUE = (1e100);
constexpr int NUM_RESOURCES = resourceNames.size();
constexpr int TOTAL_SECONDS = ((EVENT_DURATION_DAYS)*24*3600+(EVENT_DURATION_HOURS)*3600+(EVENT_DURATION_MINUTES)*60+EVENT_DURATION_SECONDS);

map<int, string> upgradeNames;
bool allInsertsWorse = false;
bool allRemovesWorse = false;
bool allSwapsWorse = false;
bool allTrimsWorse = false;
bool allRotationsWorse = false;

std::array<double, TOTAL_SECONDS> timeNeededSeconds{};

// END PROGRAM SETTINGS -----------------------------------------------------------------
// UTILITY FUNCTIONS --------------------------------------------------------------------
template <typename T>
void printVector(vector<T>& x, std::ostream& out = std::cout) {
    for (auto item : x) out << item << ",";
    out << endl;
}
class Logger {
    int interval;
    mutable int counter = 0;
public:
    Logger(int outputInterval) : interval(outputInterval) {}

    void logImprovement(const std::string& type, std::vector<int>& path, const double score) const {
        if (++counter >= interval) {
            std::cout << "Improved path " << type << " ";
            printVector(path);
            std::cout << "Score: " << score << endl;
            counter = 0;
        }
    }
};
struct SearchContext {
    Logger logger;
    const vector<double>& resources;
    const vector<int>& levels;
    mt19937& randomEngine;
    int totalSeconds;
};
void nameUpgrades() {
    for (int i = 0; i < NUM_RESOURCES; i++) {
        upgradeNames[i] = std::string(resourceNames[i]) + "_Level";
        upgradeNames[i + NUM_RESOURCES] = std::string(resourceNames[i]) + "_Speed";
    }
    upgradeNames[20] = "Complete";
}
vector <int> adjustFullPath(vector<int>& levels){
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
void printFormattedResults(vector<int>& upgradePath, vector<int>& simulationLevels, vector<double>& simulationResources, double finalScore) {
    cout << "Upgrade Path: ";
    printVector(upgradePath);
    cout << "Final Resource Counts: ";
    printVector(simulationResources);
    cout << "Final Upgrade Levels: ";
    printVector(simulationLevels);

    cout << "Event Currency: " << simulationResources[9] << endl;
    cout << "Free Exp (" << DLs << " DLs): " 
        << simulationResources[7] * (500.0 + DLs) / 5.0 
        << " (" << simulationResources[7] << " levels * cycles)" << endl;
    cout << "Pet Stones: " << simulationResources[8] << endl;
    cout << "Growth (" << UNLOCKED_PETS << " pets): " 
        << simulationResources[6] * UNLOCKED_PETS / 100.0 
        << " (" << simulationResources[6] << " levels * cycles)" << endl;
    cout << "Score: " << finalScore << endl;
}
void preprocessBusyTimes(const std::vector<double>& startHours, const std::vector<double>& endHours) {
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
vector <int> generateRandomPath(int length) {
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
        << (int)elapsedSeconds/60%60 << " minutes" << endl;
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
        1.0/1800.0, 
        1.0/2500.0, 
        1.0/1200.0, 
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
            cost[1] = baseCost;
            break;
        case 1:
            cost[1] = baseCost * 0.8;
            break;
        case 2:
            cost[1] = baseCost * 10;
            break;
        case 3:
            cost[1] = baseCost;
            cost[0] = baseCost;
            break;
        case 4:
            cost[2] = baseCost;
            break;
        case 5:
            cost[2] = baseCost;
            break;
        case 6:
            cost[0] = baseCost * 1.2;
            cost[3] = baseCost;      
            break;
        case 7:
            cost[2] = baseCost;    
            cost[4] = baseCost * 3;
            break;
        case 8:
            cost[2] = baseCost * 0.7;
            cost[5] = baseCost * 0.5;
            break;
        case 9:
            cost[3] = baseCost;
            cost[4] = baseCost;
            cost[5] = baseCost;
            break;
    }
    
    // Calculate current production rates
    productionRates[0] = levels[0] * cycleTimeMultiplier[0] * speedMultipliers[levels[10]]; // A loop is deliberately omitted here for efficiency.
    productionRates[1] = levels[1] * cycleTimeMultiplier[1] * speedMultipliers[levels[11]]; // Short loops are computationlly slower than writing out each case.
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
        if (neededResources <= 0) continue; // Already have enough resources
        
        if (productionRates[i] == 0) return INFINITY_VALUE; // Can't produce this resource yet
        
        // Find the maximum time needed across all resources to calculated expected upgrade time
        timeNeeded = max(timeNeeded, neededResources / productionRates[i]);
    }

    // Check if the user can purchase the upgrade at the specified time; given their schedule
    int busyLookupIndex = static_cast<int>(timeElapsed + timeNeeded);
    // timeNeeded += additionalTimeNeeded(timeElapsed + timeNeeded);
    if (busyLookupIndex < TOTAL_SECONDS){
        timeNeeded += timeNeededSeconds[busyLookupIndex];
    };

    // Check if we have enough time left
    if (timeNeeded >= remainingTime || upgradeType == (2 * NUM_RESOURCES)) {
        timeNeeded = remainingTime;
        resources[0] += productionRates[0] * timeNeeded; // A loop is deliberately omitted here for efficiency.
        resources[1] += productionRates[1] * timeNeeded; // Short loops are computationlly slower than writing out each case.
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

    resources[0] += productionRates[0] * timeNeeded - cost[0]; // A loop is deliberately omitted here for efficiency.
    resources[1] += productionRates[1] * timeNeeded - cost[1];
    resources[2] += productionRates[2] * timeNeeded - cost[2];
    resources[3] += productionRates[3] * timeNeeded - cost[3];
    resources[4] += productionRates[4] * timeNeeded - cost[4];
    resources[5] += productionRates[5] * timeNeeded - cost[5];
    resources[6] += productionRates[6] * timeNeeded - cost[6];
    resources[7] += productionRates[7] * timeNeeded - cost[7];
    resources[8] += productionRates[8] * timeNeeded - cost[8];
    resources[9] += productionRates[9] * timeNeeded - cost[9];

    // Perform the upgrade
    levels[upgradeType]++;
    
    return timeNeeded;
}
double simulateUpgradePath(vector<int>& path, vector<int>& levels, vector<double>& resources, double time = TOTAL_SECONDS, bool display = false) {
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
    score += resources[8] * (PET_STONES_WEIGHT);   // Pet Stones
    score += resources[6] * (GROWTH_WEIGHT);       // Growth
    
    return score;
}
double evaluatePath(vector<int>& path, const SearchContext& context){
    thread_local vector<double> testResources = context.resources;
    thread_local vector<int> testLevels = context.levels;
    testResources = context.resources;
    testLevels = context.levels;
    simulateUpgradePath(path, testLevels, testResources, context.totalSeconds);
    return calculateScore(testResources);
}
bool tryInsertUpgrade(vector<int>& upgradePath, double& finalScore, const SearchContext& context) {

    int pathLength = upgradePath.size();
    thread_local vector<int> candidatePath;

    uniform_int_distribution<> positionDist(0, pathLength);
    int startPosition = positionDist(context.randomEngine);
    const int maxTypes = (allowSpeedUpgrades ? NUM_RESOURCES * 2 : NUM_RESOURCES);

    for (int i = 0; i < pathLength; i++) {
        candidatePath = upgradePath;
        int modulatedInsertPosition = (i + startPosition) % pathLength;
        candidatePath.insert(candidatePath.begin() + modulatedInsertPosition, 0);

        int startingUpgradeType = context.randomEngine() % (allowSpeedUpgrades ? NUM_RESOURCES * 2 : NUM_RESOURCES);
        for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {

            upgradeType = (upgradeType + startingUpgradeType) % maxTypes;
            candidatePath[modulatedInsertPosition] = upgradeType;

            double testScore = evaluatePath(candidatePath, context);

            if (testScore > finalScore) {
                upgradePath.insert(upgradePath.begin() + modulatedInsertPosition, upgradeType);
                finalScore = testScore;
                context.logger.logImprovement("Insert", upgradePath, finalScore);
                return true;
            }
        }
    }
    allInsertsWorse = true;
    return false;
}
bool tryRemoveUpgrade(vector<int>& upgradePath, double& finalScore, const SearchContext& context) {

    int pathLength = upgradePath.size() - 1;
    thread_local vector<int> candidatePath;
    candidatePath = upgradePath;

    uniform_int_distribution<> swapDist(0, pathLength - 2);
    int startPos = swapDist(context.randomEngine);

    for (int i = 0; i < pathLength; i++) {
        int removePos = (i + startPos) % (pathLength);

        if (!allowSpeedUpgrades && upgradePath[removePos] >= NUM_RESOURCES) continue;

        candidatePath = upgradePath;
        candidatePath.erase(candidatePath.begin() + removePos);
        double testScore = evaluatePath(candidatePath, context);

        if (testScore >= finalScore) {
            finalScore = testScore;
            upgradePath = candidatePath;
            context.logger.logImprovement("Remove", upgradePath, finalScore);
            return true;
        }
    }
    allRemovesWorse = true;
    return false;

}
bool tryTrimPath(vector<int>& upgradePath, double& finalScore, const SearchContext& context) {

    int pathLength = upgradePath.size() - 1;
    thread_local vector<int> candidatePath;
    candidatePath = upgradePath;
    double bestNewScore = finalScore;
    int bestRemovePosition = -1;

    candidatePath.pop_back();

    for (int i = pathLength - 2; i >= 0; i--) {
        vector<int> shortenedPath(upgradePath.begin(), upgradePath.begin() + i + 1);
        shortenedPath.push_back(NUM_RESOURCES * 2);

        if (!allowSpeedUpgrades && upgradePath[i] >= NUM_RESOURCES) continue;

        double testScore = evaluatePath(candidatePath, context);

        if (testScore >= bestNewScore) {
            bestNewScore = testScore;
            bestRemovePosition = i;
        }
    }

    if (bestNewScore >= finalScore && bestRemovePosition > 0) {
        vector<int> newPath(upgradePath.begin(), upgradePath.begin() + bestRemovePosition + 1);
        newPath.push_back(NUM_RESOURCES * 2);

        upgradePath = newPath;
        finalScore = bestNewScore;
        context.logger.logImprovement("Trim", upgradePath, finalScore);
        return true;
    }
    allTrimsWorse = true;
    return false;
}
bool trySwapUpgrades(vector<int>& upgradePath, double& finalScore, const SearchContext& context) {

    int pathLength = upgradePath.size() - 1;
    thread_local vector<int> candidatePath;
    candidatePath = upgradePath;
    double newScore = finalScore;

    uniform_int_distribution<> swapDist(0, pathLength - 2);
    int startPos = swapDist(context.randomEngine);

    for (int i2 = 0; i2 < pathLength - 1; i2++) { 
        for (int j2 = i2 + 1; j2 < pathLength - 1; j2++) { 
            int i = (i2 + startPos) % (pathLength - 1);
            int j = (j2 + startPos) % (pathLength - 1);

            if (candidatePath[i] == candidatePath[j]) continue;

            swap(candidatePath[i], candidatePath[j]);

            double testScore = evaluatePath(candidatePath, context);

            if (newScore > finalScore) {
                upgradePath = candidatePath;
                finalScore = newScore;
                
                context.logger.logImprovement("Swap", upgradePath, finalScore);
                return true;
            }

            swap(candidatePath[i], candidatePath[j]);
        }
    }
    allSwapsWorse = true;
    return false;
}
bool tryRotateSubsequences(vector<int>& upgradePath, double& finalScore, const SearchContext& context) {

    int pathLength = upgradePath.size() - 1;
    thread_local vector<int> candidatePath;
    double newScore = finalScore;

    uniform_int_distribution<> rotateDist(0, pathLength - 3);
    int i = rotateDist(context.randomEngine);
    uniform_int_distribution<> rotateDist2(i+2, pathLength - 1);
    int j = rotateDist2(context.randomEngine);

    for (int k = 0; k < j-i; k++) {
        candidatePath = upgradePath;

        int offset = (k + 2) / 2;
        bool isLeft = (k % 2 == 0); 

        if(isLeft)  rotate(candidatePath.begin() + i, candidatePath.begin() + i + offset, candidatePath.begin() + j + 1); // Left rotation
        else        rotate(candidatePath.begin() + i, candidatePath.begin() + j - offset + 1, candidatePath.begin() + j + 1); // Right rotation

        double testScore = evaluatePath(candidatePath, context);

        if (newScore > finalScore) {
            upgradePath = candidatePath;
            finalScore = newScore;
            context.logger.logImprovement("Rotation", upgradePath, finalScore);
            return true;
        }
    }

    return false;
}
bool exhaustRotateSubsequences(vector<int>& upgradePath, double& finalScore, const SearchContext& context) {


    int pathLength = upgradePath.size() - 1;
    int maxIndex = pathLength - 1;
    thread_local vector<int> candidatePath;
    double newScore = finalScore;

    uniform_int_distribution<> rotateDist(0, maxIndex - 2);
    int i = rotateDist(context.randomEngine);

    for (int i2 = 0; i2 < maxIndex - 1; i2++){
        int i3 = (i + i2) % (maxIndex - 1);
        uniform_int_distribution<> rotateDist2(0, maxIndex-i3); 
        int j = rotateDist2(context.randomEngine);
        for (int j2 = 0; j2 < maxIndex - i3 - 1; j2++){ 
            int j3 = i3 + 2 + ((j + j2) % (maxIndex - i3 - 1));
            for (int k = 0; k < j3-i3; k++) {
                candidatePath = upgradePath;

                int offset = (k + 2) / 2;
                bool isLeft = (k % 2 == 0);

                if(isLeft)  rotate(candidatePath.begin() + i3, candidatePath.begin() + i3 + offset, candidatePath.begin() + j3 + 1);
                else        rotate(candidatePath.begin() + i3, candidatePath.begin() + j3 - offset + 1, candidatePath.begin() + j3 + 1);

                double testScore = evaluatePath(candidatePath, context);

                if (newScore > finalScore) {
                    upgradePath = candidatePath;
                    finalScore = newScore;
                    context.logger.logImprovement("Rotation", upgradePath, finalScore);
                    return true;
                }
            }

        }
    }
    allRotationsWorse = true;
    return false;
}
void optimizeUpgradePath(vector<int>& upgradePath, double& finalScore, const SearchContext& context, int maxIterations = 10000) {

    random_device seed;
    mt19937 randomEngine(seed());

    int iterationCount = 0;
    int noImprovementStreak = 0;

    while (noImprovementStreak < maxIterations) {
        iterationCount++;
        bool improved = false;

        int strategy = iterationCount % 100;

        if(allRotationsWorse){
            noImprovementStreak = maxIterations;
        }
        else if (allInsertsWorse && allRemovesWorse && allSwapsWorse && allTrimsWorse)
                                                  { improved = exhaustRotateSubsequences(upgradePath, finalScore, context);}
        else if (strategy < 15 && !allInsertsWorse) improved = tryInsertUpgrade(upgradePath, finalScore, context);
        else if (strategy < 30 && !allRemovesWorse) improved = tryRemoveUpgrade(upgradePath, finalScore, context);
        else if (strategy < 31 && !allTrimsWorse)   improved = tryTrimPath(upgradePath, finalScore, context);
        else if (strategy < 32)                     improved = tryRotateSubsequences(upgradePath, finalScore, context);
        else if (!allSwapsWorse)                    improved = trySwapUpgrades(upgradePath, finalScore, context);
        
        if (improved) {
            noImprovementStreak = 0;
            allInsertsWorse = false;
            allRemovesWorse = false;
            allSwapsWorse = false;
            allTrimsWorse = false;
            allRotationsWorse = false;
            continue;
        } 
        else noImprovementStreak++;
    }
}
int main() {

    nameUpgrades();
    preprocessBusyTimes(busyTimesStart, busyTimesEnd);
    
    upgradePath.reserve(500);
    if (upgradePath.empty()) {          
        upgradePath = generateRandomPath(TOTAL_SECONDS/3600);
    }

    vector<int>     simulationLevels(currentLevels), levelsCopy(currentLevels);
    vector<double>  simulationResources(resourceCounts);

    if (isFullPath) {
        levelsCopy = adjustFullPath(levelsCopy);
    }

    // Run a simulation of the upgrade path
    simulateUpgradePath(upgradePath, simulationLevels, simulationResources, TOTAL_SECONDS, false);
    double finalScore = calculateScore(simulationResources);
    printFormattedResults(upgradePath, simulationLevels, simulationResources, finalScore);

    // OPTIMIZATION ALGORITHM

    if (runOptimization) {
        srand(time(0));
        random_device seed;
        mt19937 randomEngine(seed());
        Logger logger(outputInterval);
        struct SearchContext context = {logger, resourceCounts, currentLevels, randomEngine, TOTAL_SECONDS};


        allInsertsWorse = false;
        allRemovesWorse = false;
        allSwapsWorse = false;
        allTrimsWorse = false;
        allRotationsWorse = false;

        cout << "Starting optimization..." << endl;
        optimizeUpgradePath(upgradePath, finalScore, context);
        cout << "Exhaustive attempts to improve the path have failed. This is likely your best path:" << endl;
        printVector(upgradePath);
        cout << "Final score: " << finalScore << endl;
    }
    return 0;
}   

