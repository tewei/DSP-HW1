#include <stdio.h>


void func(int array[]){
    array[1] = 9;
}

typedef struct CellPhone
{
    char hi[128];
} myphone;
int main(){
    int a[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    func(a);
    printf("%d\n", a[1]);
    myphone iphone;
    char *hello = iphone.hi;
    hello = "welcome to my phone";
    puts(hello);
    puts(iphone.hi);
    return 0;
}