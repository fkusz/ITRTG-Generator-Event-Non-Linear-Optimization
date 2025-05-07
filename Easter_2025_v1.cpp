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




//EDIT HERE START
#define DAYS 14 // Time until the event is over
#define HOURS 0
#define MINUITES 0
#define SECONDS 0
#define PETs # // number of unlocked pets from stats page
#define DLs # // top 50 DLs from stats page

vector<int> levels = { // Same order as in the game
					0,1,0, // (Level) Red, White, Blue
					0,0,0, // (Level) Green, Colorful, Yellow
					0,0,0, // (Level) Growth, Free Exp, Pet Stones
					0,           // (Level) Event Currency
					0,0,0, // (Speed) Red, White, Blue
					0,0,0, // (Speed) Green, Colorful, Yellow
					0,0,0, // (Speed) Growth, Free Exp, Pet Stones 
					0,     // (Speed) Event Currency
					0}; // dummy
vector<double> eggs = { // Red, White, etc.
					0, 500000, 0,
					0, 0, 0, 
					0/(PETs/100.0), 0/((500.0+DLs)/5.0), 0,
					0,
					};

int output_interval = 100; // If you want to run this code on your own device, set it to 1

//Advanced settings (or you can ignore below settings)
#define EVENT_CURRENCY_WEIGHT 100e-5 // Adjust to your own preferences
#define FREE_EXP_WEIGHT 6e-5
#define PET_STONES_WEIGHT 3.5e-5
#define GROWTH_WEIGHT 8e-5
vector<int> s = {6,8,10,10,13,18,0,3,10,13,0,6,10,5,13,3,0,16,0,13,0,3,16,0,8,13,6,15,1,13,2,3,16,18,0,6,17,0,14,13,12,16,7,1,3,16,1,5,12,6,9,19,9,2,4,19,13,15,12,2,16,8,6,1,3,18,12,16,18,2,6,8,4,6,12,14,3,5,0,6,15,2,6,9,18,19,9,7,12,8,3,17,6,18,12,6,8,15,7,2,5,6,19,2,3,18,0,8,6,2,18,3,17,4,14,0,14,0,19,15,4,6,8,9,2,4,0,6,18,17,7,2,8,6,14,9,3,14,8,2,17,6,7,4,19,8,5,2,17,6,3,14,0,19,0,8,0,17,6,4,3,7,14,15,2,6,17,9,3,19,8,7,19,6,2,1,8,17,6,9,8,7,7,6,9,8,20,}; // upgrade path
// the last element of this path must be 20
// paste the output by this simulator
bool full_path = true; // if the above path is the default path or other full path, it should be true
                       // if the above path is the output by this simulator, it should be false
bool speed_upgrade_add_rem = true; // if this flag is true, speed upgrade can be added or removed to/from the path
bool optimization = true;

// EDIT HERE END




#define INF (1e100)
#define N 10
#define TIME ((DAYS)*24*3600+(HOURS)*3600+(MINUITES)*60+SECONDS)
/*
0: Red level
1: White level
2: Blue level
3: Green level
4: Colorful level
5: Yellow level
6: Growth level
7: Free Exp level
8: Pet Stones level
9: Event Currency level
10: Red speed
...
19: Event Currency speed
20: complete
*/
map<int,string> mp2;

double update(vector<int> &levels, vector<double> &eggs, int x, double &time){
	//if(x>=N && levels[x]>=10) return 0;
	int y = x%N;
	if(x==N*2){
		y = -1;
	}
	int l = levels[x]+1;
	double multi = (3.0*l*l*l+1.0)*100.0;
	if(x>=N) multi *= 2.0;
	vector<double> cost(N, 0), p(N, 0);
	switch(y){
		case 0:
			cost[1] = multi;
			break;
		case 1:
			cost[1] = multi*0.8;
			break;
		case 2:
			cost[1] = multi*10;
			break;
		case 3:
			cost[1] = multi;
			cost[0] = multi;
			break;
		case 4:
			cost[2] = multi;
			break;
		case 5:
			cost[2] = multi;
			break;
		case 6:
			cost[0] = multi*1.2;
			cost[3] = multi;
			break;
		case 7:
			cost[2] = multi;
			cost[4] = multi*3;
			break;
		case 8:
			cost[2] = multi*0.7;
			cost[5] = multi*0.5;
			break;
		case 9:
			cost[3] = multi;
			cost[4] = multi;
			cost[5] = multi;
			break;
	}
	
	for(int i = 0; i < N; i++){
		p[i] = levels[i]*pow(1.25, min(levels[i+N], 10));
		switch(i){
			case 0:
				p[i] *= 1.0/3.0;
				break;
			case 1:
				break;
			case 2:
				p[i] *= 1.0/3.0;
				break;
			case 3:
				p[i] *= 1.0/3.0;
				break;
			case 4:
				p[i] *= 1.0/3.0;
				break;
			case 5:
				p[i] *= 1.0/3.0;
				break;
			case 6:
				p[i] *= 1.0/(30*60);
				break;
			case 7:
				p[i] *= 1.0/(41*60+40);
				break;
			case 8:
				p[i] *= 1.0/(20*60);
				break;
			case 9:
				p[i] *= 1.0/(1*3600+23*60+20);
				break;
                }
	}
	double t = 0;
	for(int i = 0; i < N; i++){
		double c = cost[i]-eggs[i];
		if(c<=0) continue;
		if(p[i]==0) return INF;
		t = max(t, c/p[i]);
	}
	if(t>=time||x==2*N){
		t = time;
		for(int i = 0; i < N; i++){
			eggs[i] += p[i]*t;
		}
		return t;
	}
	
	for(int i = 0; i < N; i++){
		eggs[i] += p[i]*t-cost[i];
	}
	levels[x]++;
	return t;
}

double simulate(vector<int> &path, vector<int> &levels, vector<double> &eggs, double time = TIME, bool display = false){
	//double time = 0;
	for(auto x: path){
		if(time<1e-3) return 0;
		double t = update(levels, eggs, x, time);
		time -= t;
		if(time<-1e-3) return -1;
		int tt = TIME-time;
		if(display) cout<<mp2[x]<<" "<<levels[x]<<" "<<(int)tt/24/3600<<" days, "<<(int)tt/3600%24<<" hours, "<<(int)tt/60%60<<" minutes"<<endl;
	}
	return time;
}

double score(vector<double> &eggs, bool display = false){
	double x = 0;
	for(int i = 0; i < N; i++) x += eggs[i]*1e-10;
	x += (min(eggs[9], 10000.0)/10000.0+max(0.0, (eggs[9]-10000))/10000*0.01)*10000*(EVENT_CURRENCY_WEIGHT); // event currency
	x += eggs[7]*(FREE_EXP_WEIGHT); // free exp
	x += eggs[8]*(PET_STONES_WEIGHT); // pet stone
	x += eggs[6]*(GROWTH_WEIGHT); // growth
	if(display){
		cout<<"SCORE WEIGHT DISPLAY"<<endl;
		cout<<EVENT_CURRENCY_WEIGHT<<endl;
		cout<<(500.0+DLs)/5.0/1e8<<endl;
		cout<<1.0/100000<<endl;
		cout<<PETs/100.0/10000<<endl;
	}
	return x;
}

template <typename T>
void printv(vector<T> &x){
	for(auto a: x) cout<<a<<",";
	cout<<endl;
}

int main(){
	mp2[0] = "Red_Level";
	mp2[1] = "White_Level";
	mp2[2] = "Blue_Level";
	mp2[3] = "Green_Level";
	mp2[4] = "Colorful_Level";
	mp2[5] = "Yellow_Level";
	mp2[6] = "GROWTH_Level";
	mp2[7] = "FREE_EXP_Level";
	mp2[8] = "PET_STONES_Level";
	mp2[9] = "EVENT_CURRENCY_Level";
	
	mp2[10] = "Red_Speed";
	mp2[11] = "White_Speed";
	mp2[12] = "Blue_Speed";
	mp2[13] = "Green_Speed";
	mp2[14] = "Colorful_Speed";
	mp2[15] = "Yellow_Speed";
	mp2[16] = "GROWTH_Speed";
	mp2[17] = "FREE_EXP_Speed";
	mp2[18] = "PET_STONES_Speed";
	mp2[19] = "EVENT_CURRENCY_Speed";
	
	mp2[20] = "Complete";
	
	vector<int> levels3(levels), lc(levels);
	vector<double> eggs3(eggs);

	if(full_path){
		lc[0]--;
		auto it = s.begin();
		while(it!=s.end()){
			if(lc[*it]>0){
				lc[*it]--;
				it = s.erase(it);
			} else{
				++it;
			}
		}
	}
	random_device seed_gen;
	mt19937 engine(seed_gen());
	
	double t = simulate(s, levels3, eggs3, TIME, true);
	printv(s);
	//cout<<(int)t/24/3600<<" days, "<<(int)t/3600%24<<" hours, "<<(int)t/60%60<<" minutes"<<endl;
	printv(eggs3);
	printv(levels3);
	cout<<"Event Currency: "<<eggs3[9]<<endl;
	cout<<"Free Exp ("<<DLs<<" DLs): "<<eggs3[7]*(500.0+DLs)/5.0<<" ("<<eggs3[7]<<" levels * cycles)"<<endl;
	cout<<"Pet Stones: "<<eggs3[8]<<endl;
	cout<<"Growth ("<<PETs<<" pets): "<<eggs3[6]*PETs/100.0<<" ("<<eggs3[6]<<" levels * cycles)"<<endl;
	t = score(eggs3);
	cout<<"Score: "<<t<<endl;
	if(!optimization) return 0;
	
	int count = 0, count2 = 1;
	while(true){
		count++;
		count2++;
		if(count2>output_interval) count2 = 1;
		double t2 = t;
		int n = s.size()-1;
		vector<int> s2(s);
		if(count%15==1){
			double t3 = 0;
			double k = 0;
			uniform_int_distribution<> dist(0, n);
			int pos = dist(engine);
			s2.insert(s2.begin()+pos, 0);
			for(int i = 0; i < (speed_upgrade_add_rem?N*2:N); i++){
				//cout<<i<<endl;
				//if(i==1) continue;
				s2[pos] = i;
				vector<double> eggs2(eggs);
				vector<int> levels2(levels);
				simulate(s2, levels2, eggs2,  TIME);
				t2 = score(eggs2);
				if(t2>t3){
					t3 = t2;
					k = i;
				}
			}
			s2[pos] = k;
			t2 = t3;
			if(t2>t){
				t = t2;
				s.swap(s2);
			} else {
			}
		} else if(count%15==3){
			double t3 = 0;
			double k = n-1;
			s2.pop_back();
			for(int i = n-1; i >= 0; i--){
				s2[i] = s[i+1];
				if(!speed_upgrade_add_rem&&s[i]>=N) continue;
				vector<double> eggs2(eggs);
				vector<int> levels2(levels);
				simulate(s2, levels2, eggs2,  TIME);
				t2 = score(eggs2);
				if(t2>t3){
					t3 = t2;
					k = i;
				}
			}
			for(int i = 0; i < k; i++){
				s2[i] = s[i];
			}
			t2 = t3;
			if(t2>t){
				t = t2;
				s.swap(s2);
			} else {
			}
		} else if(count%100!=0){
			uniform_int_distribution<> dist(0, n-2);
			int k = dist(engine);
			for(int i2 = 0; i2 < n-1; i2++){
				for(int j2 = i2+1; j2 < n-1; j2++){
					int i = (i2+k)%(n-1);
					int j = (j2+k)%(n-1);
					if(s2[i]==s2[j]) continue;
					swap(s2[i], s2[j]);
					vector<double> eggs2(eggs);
					vector<int> levels2(levels);
					simulate(s2, levels2, eggs2,  TIME);
					t2 = score(eggs2);
					if(t2>t) break;
					swap(s2[i], s2[j]);
				}
				if(t2>t) break;
			}
			if(t2>t){
				t = t2;
				s.swap(s2);
				if(count2==output_interval){
					printv(s);
					cout<<t<<endl;
				}
			} else {
				count = -1;
				//break;
			}
		} else {
			uniform_int_distribution<> dist(0, n-2);
			for(int i = dist(engine); i < n-1; i++){
				for(int j = i+1; j < n-1; j++){
					for(int k = 0; k < j-i; k++){
						rotate(s2.begin()+i, s2.begin()+i+1, s2.begin()+j+1);
						vector<double> eggs2(eggs);
						vector<int> levels2(levels);
						simulate(s2, levels2, eggs2,  TIME);
						t2 = score(eggs2);
						if(t2>t) break;
					}
					if(t2>t) break;
					rotate(s2.begin()+i, s2.begin()+i+1, s2.begin()+j+1);
				}
				if(t2>t) break;
			}
			if(t2>t){
				t = t2;
				s.swap(s2);
				if(true){
					printv(s);
					cout<<t<<endl;
				}
			} else {
			}
		}
	}
	printv(s);
	cout<<t<<endl;
	
	return 0;
}