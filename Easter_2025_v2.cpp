// Original Author: syokora
// v2 Update by: Frolf
    // added "Busy Times" feature so purchases won't be planned for when you're at work, school, sleeping, etc. Should help make the sims "one and done"
    // added additional methods for the optimizer to escape local optima (randomize subsequences, delete upgrades, etc.)
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
    #include <ctime> 
    #include <array>
    using namespace std;
    typedef long long ll;
    
    //------------------------------------------------------------------
    // USER CONFIGURATION - EDIT THESE VALUES TO MATCH YOUR GAME STATUS
    //------------------------------------------------------------------
    const int EVENT_DURATION_DAYS = 14;  // Time until the event is over
    const int EVENT_DURATION_HOURS = 0;
    const int EVENT_DURATION_MINUTES = 0;
    const int EVENT_DURATION_SECONDS = 0;
    const int UNLOCKED_PETS = 37;        // Number of unlocked pets from stats page
    const int DLs = 369;       // Top 50 DLs from stats page (Dungeon Levels)
    
    // Current upgrade levels (in same order as in game)
    vector<int> currentLevels = { 
        // Production Levels
        0, 1, 0,     // Red, White, Blue
        0, 0, 0,      // Green, Colorful, Yellow
        0, 0, 0,      // Growth, Free Exp, Pet Stones
        0,            // Event Currency
        
        // Speed Levels
        0, 0, 0,     // Red, White, Blue
        0, 0, 0,      // Green, Colorful, Yellow
        0, 0, 0,      // Growth, Free Exp, Pet Stones 
        0,            // Event Currency
        
        0             // Dummy placeholder
    };
    
    // Current resource counts (Red, White, Blue, etc.)
    vector<double> resourceCounts = { 
        0,  500000, 0,    // Red, White, Blue
        0,  0, 0,   // Green, Colorful, Yellow
        0/(UNLOCKED_PETS/100.0),   // Growth (scaled by pets)
        0/((500.0+DLs)/5.0), // Free Exp (scaled by Dungeon Levels)
        0,                        // Pet Stones
        0,                         // Event Currency
    };
    
    int outputInterval = 100; // How often to output results during optimization. Set this to 1 if running locally; or 100 if using an online compiler
    
    //------------------------------------------------------------------
    // REWARD PRIORITY SETTINGS (Higher weights incentivize more of that reward by the end of the event)
    //------------------------------------------------------------------
    const double EVENT_CURRENCY_WEIGHT = 100e-5; // Priority weight for Event Currency
    const double FREE_EXP_WEIGHT = 6e-5;     // Priority weight for Free EXP
    const double PET_STONES_WEIGHT = 3.5e-5;   // Priority weight for Pet Stones
    const double GROWTH_WEIGHT = 8e-5;       // Priority weight for Growth
    
    // Initial upgrade path sequence
    // Each number represents an upgrade type (0-9 for production levels, 10-19 for speed upgrades)
    // The last element must be 20 (indicates completion)
    vector<int> upgradePath = {};
    
    // Configuration flags
    bool isFullPath = true;        // If true, treat upgradePath as complete; if false, treat as output from simulator
    bool allowSpeedUpgrades = true; // If true, speed upgrades can be added/removed during optimization
    bool runOptimization = true;   // If true, run optimization algorithm to improve the path
    
    //------------------------------------------------------------------
    // END OF USER CONFIGURATION - BEGIN ADVANCED SETTINGS
    //------------------------------------------------------------------
    
    // What hours you start being unavailable to make purchases realtive to the start of the event.
    vector<double> busyTimesStart = {
        16,40,64,88,112,136,160,184,208,232,256,280,304,328
    }; // Example: I begin work 1 hour after the event starts; I start work again 24 hours later (25), I start work again another 24 hours later (49)
       // so enter 1, 25, 49, etc. Fractional hours are okay (4.25, 19.5, 37.75, 2.33333333333, etc.)
    
    //What hours you're able to make purchases again realtive to the start of the event
    vector<double> busyTimesEnd = {
        24,48,72,96,120,144,168,192,216,240,264,288,312,336 //End Times can be after the event ends.
    }; // Example: I work 8 hours, so I'm available 9 hours after the sim starts. Same for the next day, and every day.
    
    //------------------------------------------------------------------
    // END ADVANCED SETTINGS - BEGINCONSTANTS AND UTILITY DEFINITIONS
    //------------------------------------------------------------------
    
    // Define resource names (in same order as in game)
    constexpr std::array<const char*, 10> resourceNames = {"Red", "White", "Blue", "Green", "Colorful", "Yellow", "GROWTH", "FREE_EXP", "PET_STONES", "EVENT_CURRENCY"};
    
    constexpr double INFINITY_VALUE = (1e100);
    constexpr int NUM_RESOURCES = resourceNames.size();
    constexpr int TOTAL_SECONDS = ((EVENT_DURATION_DAYS)*24*3600+(EVENT_DURATION_HOURS)*3600+(EVENT_DURATION_MINUTES)*60+EVENT_DURATION_SECONDS);
    
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
    // Create flags that indicate if a search has been ehaustive and should be skipped later
    bool allInsertsWorse = false;
    bool allRemovesWorse = false;
    bool allSwapsWorse = false;
    
    /**
      * Calculates additional time for the upgrade to be purchased if the user is busy at the expected upgrade time
      * @param startTimes Start times of busy times
      * @param endTimes End times of busy times
      * @param expectedUpgradeTime Time when the upgrade is expected to be purchased
      * @return Time needed to purchase the upgrade
      */

    std::array<double, TOTAL_SECONDS> timeNeededSeconds{};
    
    // Convert busy times from hours to seconds and build a lookup table for blazing fast speeds
    void preprocessBusyTimes(const std::vector<double>& startHours, const std::vector<double>& endHours) {
        // Mark busy seconds
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
    
    // Check if the given time (in seconds) is within a busy period, and how much longer it will remain busy
    inline double additionalTimeNeeded(double expectedTimeSeconds) {
        int idx = static_cast<int>(expectedTimeSeconds);
        if (idx >= TOTAL_SECONDS) return 0.0;
        return timeNeededSeconds[idx];
    }

    /**
     * Interweave two paths together
     * @param path1 First path to interweave
     * @param path2 Second path to interweave
     * @return Interweaved path with the remianing of the longer path appended to the end
     */
    vector<int> mergePaths(vector<int>& path1, vector<int>& path2) {
        vector<int> mergedPath;
        size_t i = 0;

        // Interweave elements from both vectors until we reach the completion marker of either one
        while (i < path1.size()-1 && i < path2.size()-1) {
            mergedPath.push_back(path1[i]);
            mergedPath.push_back(path2[i]);
            i++;
        }
        
        // Append any remaining elements from path1
        while (i < path1.size()-1) {
            mergedPath.push_back(path1[i]);
            i++;
        }
        
        // Append any remaining elements from path2
        while (i < path2.size()-1) {
            mergedPath.push_back(path2[i]);
            i++;
        }
        
        // Add completion marker at the end
        mergedPath.push_back(NUM_RESOURCES * 2);
        return mergedPath;
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
        constexpr double cycleTimeMultiplier[10] = {
            1.0/3.0,    // 0
            1.0,        // 1
            1.0/3.0,    // 2
            1.0/3.0,    // 3
            1.0/3.0,    // 4
            1.0/3.0,    // 5
            1.0/1800.0, // 6
            1.0/2500.0, // 7
            1.0/1200.0, // 8
            1.0/5000.0  // 9
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
        
        // Check if this is the completion marker
        int resourceType = upgradeType % NUM_RESOURCES;
        if (upgradeType == NUM_RESOURCES * 2) {
            resourceType = -1;
        }
        
        int newLevel = levels[upgradeType] + 1;
        double baseCost = (3.0 * newLevel * newLevel * newLevel + 1.0) * 100.0; // 100(3x^3+1)
    
        if (upgradeType >= NUM_RESOURCES) baseCost *= 2.0;     // Speed upgrades cost twice as much

        
        double cost[10] = {};
        double productionRates[10] = {};
        
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
        productionRates[0] = levels[0] * cycleTimeMultiplier[0] * speedMultipliers[min(levels[10], 10)]; // A loop is deliberately omitted here for efficiency.
        productionRates[1] = levels[1] * cycleTimeMultiplier[1] * speedMultipliers[min(levels[11], 10)]; // Short loops with regular array indexing are slower than writing out each case.
        productionRates[2] = levels[2] * cycleTimeMultiplier[2] * speedMultipliers[min(levels[12], 10)];
        productionRates[3] = levels[3] * cycleTimeMultiplier[3] * speedMultipliers[min(levels[13], 10)];
        productionRates[4] = levels[4] * cycleTimeMultiplier[4] * speedMultipliers[min(levels[14], 10)];
        productionRates[5] = levels[5] * cycleTimeMultiplier[5] * speedMultipliers[min(levels[15], 10)];
        productionRates[6] = levels[6] * cycleTimeMultiplier[6] * speedMultipliers[min(levels[16], 10)];
        productionRates[7] = levels[7] * cycleTimeMultiplier[7] * speedMultipliers[min(levels[17], 10)];
        productionRates[8] = levels[8] * cycleTimeMultiplier[8] * speedMultipliers[min(levels[18], 10)];
        productionRates[9] = levels[9] * cycleTimeMultiplier[9] * speedMultipliers[min(levels[19], 10)];
        
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
        timeNeeded += additionalTimeNeeded(timeElapsed + timeNeeded);
        
        // Check if we have enough time left
        if (timeNeeded >= remainingTime || upgradeType == (2 * NUM_RESOURCES)) {
            timeNeeded = remainingTime;
            resources[0] += productionRates[0] * timeNeeded; // A loop is deliberately omitted here for efficiency.
            resources[1] += productionRates[1] * timeNeeded; // Short loops with regular array indexing is slower than writing out each case.
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
            
            // if (time < -1e-3) return -1;
            
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

    bool tryInsertUpgrade(vector<int>& upgradePath, double& finalScore, 
        const vector<double>& resourceCounts, const vector<int>& currentLevels,
        mt19937& randomEngine, bool allowSpeedUpgrades, bool verbose = false) {

        int pathLength = upgradePath.size() - 1;
        static vector<int> candidatePath;
        candidatePath = upgradePath;
        static vector<double> testResources;
        static vector<int> testLevels;

        // Choose random starting position for insertion
        uniform_int_distribution<> positionDist(0, pathLength);
        int insertPosition = positionDist(randomEngine);

        for (int i2 = 0; i2 < pathLength - 1; i2++) {
            candidatePath = upgradePath;
            int i = (i2 + insertPosition) % (pathLength - 1);
            // Try inserting each possible upgrade type
            candidatePath.insert(candidatePath.begin() + i, 0);

            // Test each upgrade type
            int maxTypes = (allowSpeedUpgrades ? NUM_RESOURCES * 2 : NUM_RESOURCES);
            for (int upgradeType = 0; upgradeType < maxTypes; upgradeType++) {
                candidatePath[i] = upgradeType;

                testResources = resourceCounts;
                testLevels = currentLevels;

                simulateUpgradePath(candidatePath, testLevels, testResources, TOTAL_SECONDS);
                double testScore = calculateScore(testResources);

                if (testScore > finalScore) {
                    upgradePath.insert(upgradePath.begin() + i, upgradeType);
                    finalScore = testScore;
                    if (verbose) {
                        cout << "Improved path (insert): ";
                        printVector(upgradePath);
                        cout << "New score: " << finalScore << endl;
                    }
                    return true;
                }
            }
        }

        if(verbose){
            cout << "No available inserts" << endl;
        }

        allInsertsWorse = true;
        return false;
    }

    // Helper function to try removing an upgrade and return if improvement was made
    bool tryRemoveUpgrade(vector<int>& upgradePath, double& finalScore, 
        const vector<double>& resourceCounts, const vector<int>& currentLevels,
        mt19937& randomEngine, bool allowSpeedUpgrades, bool verbose = false) {

        int pathLength = upgradePath.size() - 1;
        static vector<int> candidatePath;
        candidatePath = upgradePath;
        double bestNewScore = finalScore;
        int bestRemovePosition = -1;

        static vector<double> testResources;
        static vector<int> testLevels;

        // Try removing the completion marker temporarily
        candidatePath.pop_back();

        // Try moving the completion marker earlier
        for (int i = pathLength - 2; i >= 0; i--) {
            // Create a shortened path that ends at position i
            vector<int> shortenedPath(upgradePath.begin(), upgradePath.begin() + i + 1);
            shortenedPath.push_back(NUM_RESOURCES * 2); // Add completion marker

            // Skip speed upgrades if not allowed by allowSpeedUpgrades
            if (!allowSpeedUpgrades && upgradePath[i] >= NUM_RESOURCES) continue;

            testResources = resourceCounts;
            testLevels = currentLevels;

            simulateUpgradePath(shortenedPath, testLevels, testResources, TOTAL_SECONDS);
            double testScore = calculateScore(testResources);

            if (testScore >= bestNewScore) {
                bestNewScore = testScore;
                bestRemovePosition = i;
            }
        }

        // If we found a better shortened path
        if (bestNewScore >= finalScore && bestRemovePosition > 0) {
            // Create the new optimized path
            vector<int> newPath(upgradePath.begin(), upgradePath.begin() + bestRemovePosition + 1);
            newPath.push_back(NUM_RESOURCES * 2); // Add completion marker

            upgradePath = newPath;
            finalScore = bestNewScore;

            if (verbose) {
                cout << "Improved path (remove): ";
                printVector(upgradePath);
                cout << "New score: " << finalScore << endl;
            }

            return true;
        }
        if(verbose){
            cout << "No available removes" << endl;
        }
        allRemovesWorse = true;
        return false;
    }

    // Helper function to try swapping upgrades and return if improvement was made
    bool trySwapUpgrades(vector<int>& upgradePath, double& finalScore, 
        const vector<double>& resourceCounts, const vector<int>& currentLevels,
        mt19937& randomEngine, bool allowSpeedUpgrades, bool verbose = false) {

        int pathLength = upgradePath.size() - 1;
        static vector<int> candidatePath;
        candidatePath = upgradePath;
        double newScore = finalScore;

        static vector<double> testResources;
        static vector<int> testLevels;

        // Choose a random starting point for swap attempts
        uniform_int_distribution<> swapDist(0, pathLength - 2);
        int startPos = swapDist(randomEngine);

        // Try various swap combinations
        for (int i2 = 0; i2 < pathLength - 1; i2++) { 
            for (int j2 = i2 + 1; j2 < pathLength - 1; j2++) {  // Consider nearby swaps first
                int i = (i2 + startPos) % (pathLength - 1);
                int j = (j2 + startPos) % (pathLength - 1);

                // Skip if trying to swap identical upgrades
                if (candidatePath[i] == candidatePath[j]) continue;

                // Try the swap
                swap(candidatePath[i], candidatePath[j]);

                testResources = resourceCounts;
                testLevels = currentLevels;

                simulateUpgradePath(candidatePath, testLevels, testResources, TOTAL_SECONDS);
                newScore = calculateScore(testResources);

                if (newScore > finalScore) {
                    upgradePath = candidatePath;
                    finalScore = newScore;
                    
                    if (verbose) {
                        cout << "Improved path (swap): ";
                        printVector(upgradePath);
                        cout << "New score: " << finalScore << endl;
                    }
                    
                    return true;
                }

                // Revert swap if not better
                swap(candidatePath[i], candidatePath[j]);
            }
        }
        if(verbose){
            cout << "No available swaps" << endl;
        }
        allSwapsWorse = true;
        return false;
    }

    // Tries to rotate subsequences and return if improvement was made. Should be adapted to perform exhaustive rotations like the rest of the functions.
    bool tryRotateSubsequences(vector<int>& upgradePath, double& finalScore, 
        const vector<double>& resourceCounts, const vector<int>& currentLevels,
        mt19937& randomEngine, bool allowSpeedUpgrades, bool verbose = false) {

        int pathLength = upgradePath.size() - 1;
        static vector<int> candidatePath;
        double newScore = finalScore;

        static vector<double> testResources;
        static vector<int> testLevels;

        uniform_int_distribution<> rotateDist(0, pathLength - 4); // A maximum of three less than the final index of the path (Final index is end-marker, which shouldn't be rotated)

        // Determine start and end of a subsequence.
        int i = rotateDist(randomEngine);
        uniform_int_distribution<> rotateDist2(i+2, pathLength - 2); // Ensure the subsequence is AT LEAST THREE. Two would be a swap, which we already have a function for.
        int j = rotateDist2(randomEngine);

        // Perform all possible rotations of the subsequence. ~1-8 left or ~1-8 right rotations most commonly yield improvements, so we alternate between them before trying longer ones.
        for (int k = 0; k < j-i; k++) {
            candidatePath = upgradePath;

            int offset = (k + 1) / 2; // 1,1,2,2,3,3,...
            bool isLeft = (k % 2 == 1); // true, false, true, false,...

            if(isLeft)  rotate(candidatePath.begin() + i, candidatePath.begin() + i + offset, candidatePath.begin() + j); // Left rotation
            else        rotate(candidatePath.begin() + i, candidatePath.begin() + j - offset, candidatePath.begin() + j); // Right rotation

            testResources = resourceCounts;
            testLevels = currentLevels;

            simulateUpgradePath(candidatePath, testLevels, testResources, TOTAL_SECONDS);
            newScore = calculateScore(testResources);

            if (newScore > finalScore) {
                upgradePath = candidatePath;
                finalScore = newScore;

                if (verbose) {
                    cout << "Improved path (rotation): ";
                    printVector(upgradePath);
                    cout << "New score: " << finalScore << endl;
                }

                return true;
            }
        }

        // Reset rotation if no improvement
        candidatePath = upgradePath;

        return false;
    }

    bool tryNewPath(vector<int>& upgradePath, const vector<double>& resourceCounts, const vector<int>& currentLevels, double& finalScore);

    // Main optimization loop - to replace the current while(true) loop in main()
    void optimizeUpgradePath(vector<int>& upgradePath, double& finalScore, 
        const vector<double>& resourceCounts, const vector<int>& currentLevels,
        int maxIterations = 10000, int outputInterval = 100, bool skipAuxPath = true) {

        random_device seed;
        mt19937 randomEngine(seed());
        bool allowSpeedUpgrades = true; // Set based on your configuration

        int iterationCount = 0;
        int outputCount = 0;
        int noImprovementStreak = 0;

        while (noImprovementStreak < maxIterations) {
            iterationCount++;
            outputCount++;
            bool improved = false;

            // Modular approach based on iteration count - tune these percentages based on effectiveness
            int strategy = iterationCount % 100;

            if (allInsertsWorse && allRemovesWorse && allSwapsWorse) {
                improved = tryRotateSubsequences(upgradePath, finalScore, resourceCounts, currentLevels, 
                                            randomEngine, allowSpeedUpgrades, outputCount >= outputInterval);
            }
            else if (strategy < 15 && !allInsertsWorse) {
                // Insert strategy (15% of normal iterations)
                improved = tryInsertUpgrade(upgradePath, finalScore, resourceCounts, currentLevels, 
                                            randomEngine, allowSpeedUpgrades, outputCount >= outputInterval);
            } 
            else if (strategy < 30 && !allRemovesWorse) {
                // Remove strategy (15% of normal iterations)
                improved = tryRemoveUpgrade(upgradePath, finalScore, resourceCounts, currentLevels, 
                                            randomEngine, allowSpeedUpgrades, outputCount >= outputInterval);
            }
            else if (strategy < 31) {
                // Rotate strategy (1% of normal iterations)
                improved = tryRotateSubsequences(upgradePath, finalScore, resourceCounts, currentLevels, 
                                            randomEngine, allowSpeedUpgrades, outputCount >= outputInterval);
            }
            else if (!allSwapsWorse){
                // Swap strategy (63% of iterations - most common)
                improved = trySwapUpgrades(upgradePath, finalScore, resourceCounts, currentLevels, 
                                            randomEngine, allowSpeedUpgrades, outputCount >= outputInterval);
            }
            if (improved) {
                noImprovementStreak = 0;
                allInsertsWorse = false;
                allRemovesWorse = false;
                allSwapsWorse = false;

                // Reset output counter if we showed an update
                if (outputCount >= outputInterval) {
                    outputCount = 0;
                }
                continue;
            } 
            else {
                noImprovementStreak++;
            }
            // Perform extra rotations when stuck in local optimum
            if (noImprovementStreak % 10000 == 0 && !skipAuxPath) {
                cout << "Exhaustive attempts to improve the path have failed. This is likely your best path:" << endl;
                printVector(upgradePath);
                cout << "Final score: " << finalScore << endl;
                cout << "The program will now silently try to find a better path by starting over from scratch." << endl;
                cout << "If using an online compiler, note that this auxiliary search will not finish due to time constraints." << endl;
                bool escaped = tryNewPath(upgradePath, resourceCounts, currentLevels, finalScore);

                if (escaped) {
                    noImprovementStreak = 0;
                    allInsertsWorse = false;
                    allRemovesWorse = false;
                    allSwapsWorse = false;
                }
            }
        }
        

    // Final output
        if (skipAuxPath) {
            cout << "Auxiliary optimization complete" << endl;
        }
        else {
            cout << "Optimization complete after " << iterationCount << " iterations." << endl;
        }
    }
    
    bool tryNewPath(vector<int>& upgradePath, const vector<double>& resourceCounts, const vector<int>& currentLevels, double& finalScore){
        static vector<int> auxiliaryPath;
        static vector<double> testResources;
        static vector<int> testLevels;
        double auxiliaryScore = 0;

        auxiliaryPath = {};
        testResources = resourceCounts;
        testLevels = currentLevels;

        for (int i = 0; i < TOTAL_SECONDS/10800; i++) {
            auxiliaryPath.push_back(rand() % 3 + 10 * (rand() % 2));
        }
        for (int i = 0; i < TOTAL_SECONDS/10800; i++) {
            auxiliaryPath.push_back(rand() % 3 + 10 * (rand() % 2) + 3);
        }
        for (int i = 0; i < TOTAL_SECONDS/10800; i++) {
            auxiliaryPath.push_back(rand() % 4 + 10 * (rand() % 2) + 6);
        }
        auxiliaryPath.push_back(NUM_RESOURCES * 2);

        cout << "Starting auxiliary optimization..." << endl;
        optimizeUpgradePath(auxiliaryPath, auxiliaryScore, resourceCounts, currentLevels, 10000, 1000000000);

        simulateUpgradePath(auxiliaryPath, testLevels, testResources, TOTAL_SECONDS);
        auxiliaryScore = calculateScore(testResources);

        if (auxiliaryScore > finalScore) {
            upgradePath = auxiliaryPath;
            finalScore = auxiliaryScore;

            cout << "Improved path (Auxiliary): ";
            printVector(upgradePath);
            cout << "New score: " << finalScore << endl;

            return true;
        }
        else {
            auxiliaryPath = mergePaths(upgradePath, auxiliaryPath);
            optimizeUpgradePath(auxiliaryPath, auxiliaryScore, resourceCounts, currentLevels, 10000, 1000000000);
            testResources = resourceCounts;
            testLevels = currentLevels;

            simulateUpgradePath(auxiliaryPath, testLevels, testResources, TOTAL_SECONDS);
            auxiliaryScore = calculateScore(testResources);

            if (auxiliaryScore > finalScore) {

                upgradePath = auxiliaryPath;
                finalScore = auxiliaryScore;
                cout << "Improved path (Auxiliary Merge): ";
                printVector(upgradePath);
                cout << "New score: " << finalScore << endl;

                return true;
            }
        }
        return false;
    }
    int main() {
        srand(time(0));

        upgradePath.reserve(500);

        // Initialize upgrade names
        for (int i = 0; i < NUM_RESOURCES; i++) {
            upgradeNames[i] = std::string(resourceNames[i]) + "_Level";
            upgradeNames[i + NUM_RESOURCES] = std::string(resourceNames[i]) + "_Speed";
        }
        upgradeNames[20] = "Complete";
    
        // Convert busy times to seconds
        preprocessBusyTimes(busyTimesStart, busyTimesEnd);
    
        // If no upgrade path is provided, generate a random inital path with one upgrade per hour
        if (upgradePath.empty()) {          
            for (int i = 0; i < TOTAL_SECONDS/10800; i++) {
                upgradePath.push_back(rand() % 3 + 10 * (rand() % 2));
            }
            for (int i = 0; i < TOTAL_SECONDS/10800; i++) {
                upgradePath.push_back(rand() % 3 + 10 * (rand() % 2) + 3);
            }
            for (int i = 0; i < TOTAL_SECONDS/10800; i++) {
                upgradePath.push_back(rand() % 4 + 10 * (rand() % 2) + 6);
            }
            upgradePath.push_back(NUM_RESOURCES * 2);
        }
    
        // Create working copies of levels and resources
        vector<int> simulationLevels(currentLevels), levelsCopy(currentLevels);
        vector<double> simulationResources(resourceCounts);
    
        // Remove already bought upgrades from the upgrade path if a pruned path isn't provided from the start.
        if (isFullPath) {
            levelsCopy[1]--; // Adjust first level
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

        // Remove speed upgrades with already maxed Levels that may be leftover in the path
        auto it = upgradePath.begin();
        while (it != upgradePath.end()) {
            if (*it >= NUM_RESOURCES && levelsCopy[*it] == 10) {
                it = upgradePath.erase(it);
            } else {
                ++it;
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
    
        if (runOptimization) {
            cout << "Starting optimization..." << endl;
            optimizeUpgradePath(upgradePath, finalScore, resourceCounts, currentLevels, 10000, outputInterval, false);
        }
    
        // Final output
        cout << "Optimized upgrade path: ";
        printVector(upgradePath);
        cout << "Final score: " << finalScore << endl;
    
        return 0;
    }
    