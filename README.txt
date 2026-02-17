Using the Script

---------- Running via online compiler ---------

0) Copy the script into either a text editor or into the online compiler, make sure you update the values in the #USER SETTINGS area. Most notably, use your number of unlocked pets, and number of DLs per your stats page. Make sure to updated these as needed throughout the event. Also ensure your time until event end is adjusted as needed throughout. You may also want to adjust ADVANCED SETTINGS (elaborated below) but these are not required for a decent resilt. 
Once the values are updated, you can proceed with running the code. Ensure you are using an online compiler that lets you use C++ 23 or newer, as old compilers may run into errors. 

1) The first run will time itself out, you can either wait for this to happen or wait until the score does not change drastically. Each set of data being output is the upgrade path on top, with the respective score underneath [upgrade path is the long list of numbers separated by commas i.e. CSV].

2) If the first run is done and you want to optimize further, you can copy the comma separated values into the upgradePath on line 59 of the script. This will continue optimizing where you left off, ensuring you can optimize as long as you wish. 

3) Given enough time, the sim will end on its own (after 10,000 consecutive failed improvements). In this case, the code will output your upgrade plan in plain, readable text telling you which upgrades to buy in which order. If you don't want to wait for the sim to end on its own, you should copy an upgrade path your happy with into like 59 of the script and set line 62 to "false" and run the code again. This will skip the optimization process and simply return your upgrade plan in plain text. 

4) Another thing to note, is that the time values for each step are time from script-start. If you see 28min and the next one is 38min that means you need to wait 10minutes between those two steps.

You can use https://www.programiz.com/cpp-programming/online-compiler/ (or any other C++ compiler)

---------- Running Locally ---------

1) If running locally, ensure you are using an up to date compiler. This script was developed and tested using C++23, and it's been noted that older compiler standards will often error. The sim has an option on line 63 to run endlessly, in case you want to let a computer spend many hours finding the best possible path. Otherwise, the script will end when there have been 10,000 consecutive failures to improve the path. 

2) If you don't wait for the sim to end naturally, you need to take one of those sims with a score you're happy with and paste it back into upgradePath on line 59, recompile, and run again. This time, the results will print immediately. 

----------- Advanced Settings ------------ 

Output interval: one of the slowest part of this program is actually printing out results. You can change how often this occurs by adjusting this setting. 1000ms = 1s. I recommend a value of 3-5 seconds. This has no bearing on the accuracy of the sim; just reduces output clutter. Very small values will make the sim sluggish. 

Weights: if you want to prioritize certain event resources over others, you can change how heavily certain outputs are scored here. Make sure you keep event currency decently higher than other weights, or the sim will ignore it. Higher weights (relative to others) makes sure the sim will produce more of that resource. 

Busy Times: to account for work meetings, sleeping, and whatever else you may have happening in life, you can pre-schedule times where you are unavailable to purchase upgrades. This allows you to have the sim optimize AROUND your schedule, rather than in ignorance to it. Like the buy times for the output, the lists you create should be relative to the start of your simulation. An example below should help clarify the use: 

Let's say I am stimming to start my event at noon and plan to sleep from 10 PM to 6 AM (22:00-6:00):

busyTimesStart = {10,34,58,82...}
busyTimesEnd = {18,42,66,90...}

Telling the optimizer I'm going to be away starting in 10 hours, and I'll be available again in 18 hours. I add 24 to each value to tell it this will happen every day.



