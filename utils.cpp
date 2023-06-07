#include "utils.h"


void print_vector(vector<long long> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}

void print_vector(vector<Z_NR<mpz_t>> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}



void print_vector(vector<int> v,int index_start, int index_end){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}

void print_map(map<long long, long long >  mp){
    cout<<"{ ";
    for (map<long long, long long >  ::iterator it = mp.begin();
		it != mp.end(); it++) {
        cout<<it-> first<<":"<<it->second<<", ";
	}
    cout<<"}"<<endl;
}