// Original Author: syokora
// v2 Update by: Frolf
    // added "Busy Times" feature so purchases won't be planned for when you're at work, school, sleeping, etc. Should help make the sims "one and done"
    // (WIP) added additional methods for the optimizer to escape local optima (randomize subsequences, delete upgrades, etc.)
    // (WIP) added the option for absolute timestamps for upgrade buy-times rather than timestamps relative to the start of the sim
    // updated variable names for imporved maintainability; made some generalization changes to make it easier to update for different events;
    // cleaned up some extremely minor inefficiencies; added some output clarity

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <random>
#include <iterator>
#include <numeric>
#include <algorithm>
using namespace std;
typedef long long ll;

//------------------------------------------------------------------
// USER CONFIGURATION - EDIT THESE VALUES TO MATCH YOUR GAME STATUS
//------------------------------------------------------------------
#define EVENT_DURATION_DAYS 4  // Time until the event is over
#define EVENT_DURATION_HOURS 7
#define EVENT_DURATION_MINUTES 0
#define EVENT_DURATION_SECONDS 0
#define UNLOCKED_PETS 32        // Number of unlocked pets from stats page
#define DLs 303       // Top 50 DLs from stats page (Dungeon Levels)

// Current upgrade levels (in same order as in game)
vector<int> currentLevels = { 
    // Production Levels
    25, 18, 16,     // Red, White, Blue
    18, 5, 12,      // Green, Colorful, Yellow
    23, 0, 24,      // Growth, Free Exp, Pet Stones
    9,            // Event Currency
    
    // Speed Levels
    10, 10, 10,     // Red, White, Blue
    10, 4, 10,      // Green, Colorful, Yellow
    10, 0, 10,      // Growth, Free Exp, Pet Stones 
    9,            // Event Currency
    
    0             // Dummy placeholder
};

// Current resource counts (Red, White, Blue, etc.)
vector<double> resourceCounts = { 
    1073000,  139520, 1481000,    // Red, White, Blue
    1807000,  106776, 1042000,   // Green, Colorful, Yellow
    19603/(UNLOCKED_PETS/100.0),   // Growth (scaled by pets)
    0/((500.0+DLs)/5.0), // Free Exp (scaled by Dungeon Levels)
    85062,                        // Pet Stones
    3609,                         // Event Currency
};

int outputInterval = 1; // How often to output results during optimization. Set this to 1 if running locally; or 100 if using an online compiler

//------------------------------------------------------------------
// REWARD PRIORITY SETTINGS (Higher weights incentivize more of that reward by the end of the event)
//------------------------------------------------------------------
#define EVENT_CURRENCY_WEIGHT 100e-5 // Priority weight for Event Currency
#define FREE_EXP_WEIGHT 1e-7     // Priority weight for Free EXP
#define PET_STONES_WEIGHT 6e-5   // Priority weight for Pet Stones
#define GROWTH_WEIGHT 6e-5       // Priority weight for Growth

// Initial upgrade path sequence
// Each number represents an upgrade type (0-9 for production levels, 10-19 for speed upgrades)
// The last element must be 20 (indicates completion)
vector<int> upgradePath = {14,5,19,9,3,7,0,0,8,6,9,3,8,17,7,17,7,17,4,14,4,14,4,14,7,4,2,14,6,14,9,4,17,7,19,2,9,7,3,17,7,17,0,8,0,6,9,7,3,8,4,4,6,9,3,9,8,2,2,10,10,14,12,2,8,4,6,14,4,2,12,18,16,14,0,14,12,2,16,18,4,14,18,10,8,5,6,9,4,6,6,8,4,2,6,16,10,2,10,4,4,7,14,16,18,6,2,7,18,16,14,12,17,4,4,10,16,12,0,4,8,16,12,10,0,18,2,0,0,16,2,6,16,18,0,8,4,8,12,12,12,12,2,16,8,2,14,12,6,18,16,16,4,18,2,};

// Configuration flags
bool isFullPath = false;        // If true, treat upgradePath as complete; if false, treat as output from simulator
bool allowSpeedUpgrades = true; // If true, speed upgrades can be added/removed during optimization
bool runOptimization = true;   // If true, run optimization algorithm to improve the path

//------------------------------------------------------------------
// END OF USER CONFIGURATION - BEGIN ADVANCED SETTINGS
//------------------------------------------------------------------

// What hours you start being unavailable to make purchases realtive to the start of the event.
vector<double> busyTimesStart = {
    2,26,50,74,98,122,146,170,194,218,242,266
}; // Example: I begin work 1 hour after the event starts; I start work again 24 hours later (25), I start work again another 24 hours later (49)
   // so enter 1, 25, 49, etc. Fractional hours are okay (4.25, 19.5, 37.75, 2.33333333333, etc.)

//What hours you're able to make purchases again realtive to the start of the event
vector<double> busyTimesEnd = {
    10,34,58,82,106,130,154,178,202,226,250,274 //End Times can be after the event ends.
}; // Example: I work 8 hours, so I'm available 9 hours after the sim starts. Same for the next day, and every day.

//------------------------------------------------------------------
// END ADVANCED SETTINGS - BEGINCONSTANTS AND UTILITY DEFINITIONS
//------------------------------------------------------------------

// Define resource names (in same order as in game)
vector<string> resourceNames = {"Red", "White", "Blue", "Green", "Colorful", "Yellow", "GROWTH", "FREE_EXP", "PET_STONES", "EVENT_CURRENCY"};

#define INFINITY_VALUE (1e100)
#define NUM_RESOURCES resourceNames.size()
#define TOTAL_SECONDS ((EVENT_DURATION_DAYS)*24*3600+(EVENT_DURATION_HOURS)*3600+(EVENT_DURATION_MINUTES)*60+EVENT_DURATION_SECONDS)

/*
Resource indices:
0: Red production level
1: White production level
2: Blue production level
3: Green production level
4: Colorful production level
5: Yellow production level
6: Growth production level
7: Free Exp production level
8: Pet Stones production level
9: Event baskets production level (Event Currency)
10-19: Speed upgrades for the same resources
20: Completion marker
*/

// Initiate map to store readable names for each upgrade type. This is populated in main()
map<int, string> upgradeNames;

/**
  * Calculates additional time for the upgrade to be purchased if the user is busy at the expected upgrade time
  * @param startTimes Start times of busy times
  * @param endTimes End times of busy times
  * @param expectedUpgradeTime Time when the upgrade is expected to be purchased
  * @return Time needed to purchase the upgrade
  */
double additionalTimeNeeded(vector<double> startTimes, vector<double> endTimes, double expectedUpgradeTime) {
    double timeNeeded = 0;
    // If the expected upgrade time is during a busy time, calculate how long until the user can get to their device
    for (int i = 0; i < startTimes.size(); i++) {
        if (startTimes[i] <= expectedUpgradeTime && endTimes[i] >= expectedUpgradeTime) {
            timeNeeded = endTimes[i] - expectedUpgradeTime;
            break;
        }
    }
    return timeNeeded;
}

/**
 * Update game state after an upgrade
 * @param levels Current upgrade levels
 * @param resources Current resource counts
 * @param upgradeType Type of upgrade to perform
 * @param remainingTime Remaining time in the event
 * @return Time taken for the upgrade
 */
double performUpgrade(vector<int>& levels, vector<double>& resources, int upgradeType, double& remainingTime) {
    // Check if this is the completion marker
    int resourceType = upgradeType % NUM_RESOURCES;
    if (upgradeType == NUM_RESOURCES * 2) {
        resourceType = -1;
    }
    
    int newLevel = levels[upgradeType] + 1;
    double baseCost = (3.0 * newLevel * newLevel * newLevel + 1.0) * 100.0; // 100(3x^3+1)

    if (upgradeType >= NUM_RESOURCES) baseCost *= 2.0;     // Speed upgrades cost twice as much
    vector<double> cost(NUM_RESOURCES, 0), productionRates(NUM_RESOURCES, 0);   // Calculate costs for the upgrade
    
    // Set upgrade costs based on resource type. Some resources' costs are scaled by a constant factor of the base cost.
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
    for (int i = 0; i < NUM_RESOURCES; i++) {
        // Base production = level * speed multiplier
        productionRates[i] = levels[i] * pow(1.25, min(levels[i + NUM_RESOURCES], 10)); //1.25 comes from the inverse of 0.8 (20% speed decrease is a 25% production increase)
        
        // Apply resource-specific rate adjustments
        switch (i) {
            case 0:
            case 2: 
            case 3: 
            case 4: 
            case 5: 
                productionRates[i] *= 1.0/3.0; // 3-second base cycle
                break;
            case 1: // 1-second base cycle (no adjustment needed)
                break;
            case 6:
                productionRates[i] *= 1.0/(30*60); // 30-minute cycle
                break;
            case 7:
                productionRates[i] *= 1.0/(41*60+40); // 41:40 minute cycle (2500 seconds)
                break;
            case 8:
                productionRates[i] *= 1.0/(20*60); // 20-minute cycle
                break;
            case 9:
                productionRates[i] *= 1.0/(1*3600+23*60+20); // 1:23:20 hour cycle (5000 seconds)
                break;
        }
    }
    
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
    timeNeeded += additionalTimeNeeded(busyTimesStart, busyTimesEnd, timeElapsed + timeNeeded);
    
    // Check if we have enough time left
    if (timeNeeded >= remainingTime || upgradeType == 2 * NUM_RESOURCES) {
        timeNeeded = remainingTime;
        // Just accumulate resources with the remaining time
        for (int i = 0; i < NUM_RESOURCES; i++) {
            resources[i] += productionRates[i] * timeNeeded;
        }
        return timeNeeded;
    }
    
    // Update resource counts after waiting for the upgrade
    for (int i = 0; i < NUM_RESOURCES; i++) {
        resources[i] += productionRates[i] * timeNeeded - cost[i];
    }
    
    // Perform the upgrade
    levels[upgradeType]++;
    
    return timeNeeded;
}

/**
 * Simulate the entire upgrade path
 * @param path Sequence of upgrades to perform
 * @param levels Current upgrade levels
 * @param resources Current resource counts
 * @param time Total time available
 * @param display Whether to display detailed progress
 * @return Remaining time after all upgrades
 */
double simulateUpgradePath(vector<int>& path, vector<int>& levels, vector<double>& resources, 
                          double time = TOTAL_SECONDS, bool display = false) {
    for (auto upgradeType : path) {
        if (time < 1e-3) return 0; // No time left
        
        double timeTaken = performUpgrade(levels, resources, upgradeType, time);
        time -= timeTaken;
        
        if (time < -1e-3) return -1; // Error in time calculation
        
        // Display progress if requested (debugging option)
        if (display) {
            int elapsedSeconds = TOTAL_SECONDS - time;
            cout << upgradeNames[upgradeType] << " " << levels[upgradeType] 
                 << " " << (int)elapsedSeconds/24/3600 << " days, " 
                 << (int)elapsedSeconds/3600%24 << " hours, " 
                 << (int)elapsedSeconds/60%60 << " minutes" << endl;
        }
    }
    return time;
}

/**
 * Calculate the score for the current state
 * @param resources Current resource counts
 * @param display Whether to display score components
 * @return Overall score based on weighted rewards
 */
double calculateScore(vector<double>& resources, bool display = false) {
    double score = 0;
    
    // Base score is sum of all resources (very small weight)
    for (int i = 0; i < NUM_RESOURCES; i++) {
        score += resources[i] * 1e-15;
    }
    
    // Add weighted scores for each final reward
    // Event Curreny (capped at 10000, with diminishing returns after to incentivize a slight overshoot)
    score += (min(resources[9], 10000.0) + max(0.0, (resources[9] - 10000)) * 0.01) * (EVENT_CURRENCY_WEIGHT);
    
    // Other rewards
    score += resources[7] * (FREE_EXP_WEIGHT);     // Free EXP
    score += resources[8] * (PET_STONES_WEIGHT);   // Pet Stones
    score += resources[6] * (GROWTH_WEIGHT);       // Growth
    
    // Display scoring weights if requested
    if (display) {
        cout << "SCORE WEIGHT DISPLAY" << endl;
        cout << "Event Currency Weight: " << EVENT_CURRENCY_WEIGHT << endl;
        cout << "Free EXP Scaling: " << (500.0 + DLs) / 5.0 / 1e8 << endl;
        cout << "Pet Stones Base: " << 1.0 / 100000 << endl;
        cout << "Growth Scaling: " << UNLOCKED_PETS / 100.0 / 10000 << endl;
    }
    
    return score;
}

/**
 * Print a vector to console
 */
template <typename T>
void printVector(vector<T>& x) {
    for (auto item : x) cout << item << ",";
    cout << endl;
}

int main() {
    // Initialize upgrade names
    for (int i = 0; i < NUM_RESOURCES; i++) {
        upgradeNames[i] = resourceNames[i] + "_Level";
        upgradeNames[i + NUM_RESOURCES] = resourceNames[i] + "_Speed";
    }
    upgradeNames[20] = "Complete";

    // Convert busy times to seconds
    for (int i = 0; i < busyTimesStart.size(); i++) {
        busyTimesStart[i] *= 3600;
        busyTimesEnd[i] *= 3600;
    }

    // If no upgrade path is provides, generate a random path of 100 upgrades
    if (upgradePath.empty()) {
        for (int i = 0; i < 100; i++) {
            upgradePath.push_back(rand() % NUM_RESOURCES * 2);
        }
    }

    // Create working copies of levels and resources
    vector<int> simulationLevels(currentLevels), levelsCopy(currentLevels);
    vector<double> simulationResources(resourceCounts);

    // Remove already bought upgrades from the upgrade path if a pruned path isn't provided from the start.
    if (isFullPath) {
        levelsCopy[0]--; // Adjust first level
        auto it = upgradePath.begin();
        while (it != upgradePath.end()) {
            if (levelsCopy[*it] > 0) {
                // Remove already completed upgrades
                levelsCopy[*it]--;
                it = upgradePath.erase(it);
            } else {
                ++it;
            }
        }
    }

    // Initialize random number generator
    random_device seed;
    mt19937 randomEngine(seed());

    // Run initial simulation with current upgrade path
    double remainingTime = simulateUpgradePath(upgradePath, simulationLevels, simulationResources, TOTAL_SECONDS, true);

    // Print results
    cout << "Upgrade Path: ";
    printVector(upgradePath);

    cout << "Final Resource Counts: ";
    printVector(simulationResources);

    cout << "Final Upgrade Levels: ";
    printVector(simulationLevels);

    // Display final rewards
    cout << "Event Currency: " << simulationResources[9] << endl;
    cout << "Free Exp (" << DLs << " DLs): " 
         << simulationResources[7] * (500.0 + DLs) / 5.0 
         << " (" << simulationResources[7] << " levels * cycles)" << endl;
    cout << "Pet Stones: " << simulationResources[8] << endl;
    cout << "Growth (" << UNLOCKED_PETS << " pets): " 
         << simulationResources[6] * UNLOCKED_PETS / 100.0 
         << " (" << simulationResources[6] << " levels * cycles)" << endl;

    // Calculate overall score
    double finalScore = calculateScore(simulationResources);
    cout << "Score: " << finalScore << endl;

    // Stop here if optimization is disabled
    if (!runOptimization) return 0;

    //------------------------------------------------------------------
    // OPTIMIZATION ALGORITHM
    //------------------------------------------------------------------
    int iterationCount = 0, outputCount = 1;

    while (true) {
        iterationCount++;
        outputCount++;
        if (outputCount > outputInterval) outputCount = 1;

        double newScore = finalScore;
        int pathLength = upgradePath.size() - 1;
        vector<int> candidatePath(upgradePath);

        // Strategy 1: Insert a new upgrade (every 15 iterations)
        if (iterationCount % 15 == 1) {
            double bestNewScore = 0;
            int bestUpgradeType = 0;

            // Choose random position for insertion
            uniform_int_distribution<> positionDist(0, pathLength);
            int insertPosition = positionDist(randomEngine);

            // Try inserting each possible upgrade type
            candidatePath.insert(candidatePath.begin() + insertPosition, 0);

            // Test each upgrade type
            int maxTypes = (allowSpeedUpgrades ? NUM_RESOURCES * 2 : NUM_RESOURCES);
            for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {
                candidatePath[insertPosition] = upgradeType;

                vector<double> testResources(resourceCounts);
                vector<int> testLevels(currentLevels);

                simulateUpgradePath(candidatePath, testLevels, testResources, TOTAL_SECONDS);
                double testScore = calculateScore(testResources);

                if (testScore > bestNewScore) {
                    bestNewScore = testScore;
                    bestUpgradeType = upgradeType;
                }
            }

            // Use the best upgrade found
            candidatePath[insertPosition] = bestUpgradeType;
            newScore = bestNewScore;

            if (newScore > finalScore) {
                finalScore = newScore;
                upgradePath.swap(candidatePath);
            }
        }
        // Strategy 2: Remove an upgrade (every 15 iterations)
        else if (iterationCount % 15 == 3) {
            double bestNewScore = 0;
            int bestRemovePosition = pathLength - 1;

            // Try removing the completion marker temporarily
            candidatePath.pop_back();

            for (int i = pathLength - 1; i >= 0; i--) {
                candidatePath[i] = upgradePath[i + 1];

                // Skip speed upgrades if not allowed by allowSpeedUpgrades
                if (!allowSpeedUpgrades && upgradePath[i] >= NUM_RESOURCES) continue;

                vector<double> testResources(resourceCounts);
                vector<int> testLevels(currentLevels);

                simulateUpgradePath(candidatePath, testLevels, testResources, TOTAL_SECONDS);
                double testScore = calculateScore(testResources);

                if (testScore > bestNewScore) {
                    bestNewScore = testScore;
                    bestRemovePosition = i;
                }
            }

            // Keep original upgrades up to the best position
            for (int i = 0; i < bestRemovePosition; i++) {
                candidatePath[i] = upgradePath[i];
            }

            newScore = bestNewScore;

            if (newScore > finalScore) {
                finalScore = newScore;
                upgradePath.swap(candidatePath);
            }
        }
        // Strategy 3: Swap upgrades (most iterations)
        else if (iterationCount % 75 != 0) {
            // Choose a random starting point for swap attempts
            uniform_int_distribution<> swapDist(0, pathLength - 2);
            int startPos = swapDist(randomEngine);

            // Try various swap combinations
            for (int i2 = 0; i2 < pathLength - 1; i2++) {
                for (int j2 = i2 + 1; j2 < pathLength - 1; j2++) {
                    int i = (i2 + startPos) % (pathLength - 1);
                    int j = (j2 + startPos) % (pathLength - 1);

                    // Skip if trying to swap identical upgrades
                    if (candidatePath[i] == candidatePath[j]) continue;

                    // Try the swap
                    swap(candidatePath[i], candidatePath[j]);

                    vector<double> testResources(resourceCounts);
                    vector<int> testLevels(currentLevels);

                    simulateUpgradePath(candidatePath, testLevels, testResources, TOTAL_SECONDS);
                    newScore = calculateScore(testResources);

                    if (newScore > finalScore) break;

                    // Revert swap if not better
                    swap(candidatePath[i], candidatePath[j]);
                }
                if (newScore > finalScore) break;
            }

            if (newScore > finalScore) {
                finalScore = newScore;
                upgradePath.swap(candidatePath);

                // Output progress at specified intervals
                if (outputCount == outputInterval) {
                    cout << "Improved path: ";
                    printVector(upgradePath);
                    cout << "New score: " << finalScore << endl;
                }
            } else {
                // Reset if we're not finding improvements
                iterationCount = -1;
            }
        }
        // Strategy 4: Rotate subsequences (every 100 iterations)
        else {
            uniform_int_distribution<> rotateDist(0, pathLength - 2);

            // Try rotating various subsequences
            for (int i = rotateDist(randomEngine); i < pathLength - 1; i++) {
                for (int j = i + 1; j < pathLength - 1; j++) {
                    for (int k = 0; k < j - i; k++) {
                        // Rotate the subsequence
                        rotate(candidatePath.begin() + i, candidatePath.begin() + i + 1, candidatePath.begin() + j + 1);

                        vector<double> testResources(resourceCounts);
                        vector<int> testLevels(currentLevels);

                        simulateUpgradePath(candidatePath, testLevels, testResources, TOTAL_SECONDS);
                        newScore = calculateScore(testResources);

                        if (newScore > finalScore) break;
                    }
                    if (newScore > finalScore) break;

                    // Reset rotation if no improvement
                    rotate(candidatePath.begin() + i, candidatePath.begin() + i + 1, candidatePath.begin() + j + 1);
                }
                if (newScore > finalScore) break;
            }

            if (newScore > finalScore) {
                finalScore = newScore;
                upgradePath.swap(candidatePath);

                // Always output rotation improvements; these are crucial to break out of local optima
                cout << "Improved path (rotation): ";
                printVector(upgradePath);
                cout << "New score: " << finalScore << endl;
            }
        }
    }

    // Final output
    cout << "Optimized upgrade path: ";
    printVector(upgradePath);
    cout << "Final score: " << finalScore << endl;

    return 0;
}
