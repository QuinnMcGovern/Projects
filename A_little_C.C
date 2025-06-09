//This file is used to replace some basic Python files in C based on what I did in JCP265 and PHY426

//First I will do a simulation to find the Catalan Numbers!!!

#include <stdio.h>//This is used to import a few needed functions

int Cn = 1; //sets our inital value for the first of the numbers
int n = 0;//sets up out countre we will use to calculate higher order terms

int Results_Array[50] = {}; //This array will hold our results when we have them, we prealocate for 50 ints

//we can also come up with a function to calculate factorials
int factorial (int l) {
    /*This function is used to calculate the given value of a factorial of an inputed natural number l
    such that factorials are calculated as by output = l * (l-1) * (l-2) * (l-3) * ... ( l - l - 1). Only accurate for 
    natural numbers and zero, negative numbers are not usefull, it will just return 1.
    Input: int
    Returns int*/
    while (l >= 1) { //i our value is 1 or more we use our recurion algorithum
        return l * factorial (l - 1);
    }
    while (l <= 0) { //if our value of i is 0 or less we set the funtion to output one
        return 1;
    }
}
        

int main () {

    do
    {
        Results_Array[n] = Cn;// Saves the results from the prior loop
        Cn = ((((4.0 * n) + 2.0) / (n + 2.0))) * Cn; //This preforms the calculation to find the next value of Cn in the itteration
        n = ++n; //adds one to our countre such that we can continue to iterate

    } while (Cn < 1e9);

    int i; //defines an interger i we will use as a countre

    for (i = 0; i < 19; i++){ //we loop through the first 20 values in the array  
        printf ("%d", Results_Array[i]); // we print each of the first 20 values in the array
        printf ("\n"); //then print a blank line
    }  //confirmed to work
    
    printf ("\n"); //blanks to space out the results
    printf ("\n");

    printf ("%d",factorial(9)); //checked and works 


    return 0;
}



