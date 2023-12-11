#include<iostream>
#include<algorithm>
#include<set>
using namespace std;

void bubbleSort(int *arr, int size)
{
    int temp;
    for(int i = 0;i < size; ++ i)
    {
        for(int q = 1;q < size - i; ++ q)
        {
            if(arr[q] < arr[q - 1])
            {
                temp = arr[q];
                arr[q] = arr[q - 1];
                arr[q - 1] = temp;
            }
        }
    }
}


int main()
{
    int arr[] = {3, 4, 2, 5, 1};
    string arrr;
    arrr += "10";
    arrr += "94";
    cout << arrr << endl;
    for(int i : arr) cout << i << ' ';
    cout << endl;
    bubbleSort(arr, 5);
    for(int i : arr) cout << i << ' ';

}