#include <iostream>
using namespace :: std;
class matrix{
	public:
		matrix(){cout << "matrix init" << endl;}
        void hoge(){cout << "virtual test" << endl;}
};

class inherit : public matrix{
	int a;
	public:
	void hoge(){cout << "virtual OK" << endl; }
};

int main(){
	matrix a;
        inherit b;
	matrix::hoge();
	b::hoge();
}

