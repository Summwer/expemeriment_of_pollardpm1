#include "utils.h"



void print_map(map<long long, long long >  mp){
    cout<<"{ ";
    for (map<long long, long long >  ::iterator it = mp.begin();
		it != mp.end(); it++) {
        cout<<it-> first<<":"<<it->second<<", ";
	}
    cout<<"}"<<endl;
}