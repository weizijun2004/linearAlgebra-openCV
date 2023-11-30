#include<iostream>
#include<vector>
using namespace std;

class animal
{
public : 
    void call(string s)
    {
        cout << s << endl;
    }
    void eat()
    {
        cout << "eating" << endl;
    }
private : 
    
    string name;
    string id;
    int sex;
protected : 
    int age = 0;
};

class dog : public animal
{

public : 
    int outputAge()
    {
        cout << age << endl;
    }

private : 
    string color;
};

class dog123 : public dog
{
public : 
    void outputAge()
    {
        cout << age << endl;
    }
};

int main()
{
    vector<dog> arr1;
    dog test;
    test.call("www");
    test.eat();
    cout << test.age << endl;
    cout << "hello world." << endl;

}