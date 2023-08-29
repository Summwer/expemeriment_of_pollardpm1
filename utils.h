
#include <vector>
#include <iostream>
#include <fplll.h>

using namespace std;
using namespace fplll;


void print_map(map<long long, long long >  mp);

template <class T> inline void print_vector(vector<T> v,int index_start=0, int index_end =-1){
    cout<< "[";
    if(index_end == -1)
        index_end = v.size();
    for(int i = index_start; i < index_end; i++) {
        cout << v[i] << " ";
    }
    cout<<"]"<<endl;
}
