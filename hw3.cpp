#include <iostream>
#include <stdlib.h>
#include <cmath> // should fix my error (stupid abs couldnt handle doubles(tho didnt throw error), nice job c)
#include <omp.h>
#include <stdio.h>
#include <cstdlib>





int main(int argc, char* argv[])
{
    bool verbose = false;
    bool debug = false;
    bool showL21 = false;
    srand48(time(0));
    
    if(debug)
    {
        std::cout << "*** Debug mode ONLINE\n";
    }
    
    std::cout << "Let us do some arbitrary calculations!\n";
    
    std::cout << "We will use a matrix of size n = " << argv[1] << "\n";
    int n = atoi(argv[1]);
    
    std::cout << "We will use a threads size of p = " << argv[2] << " threads\n";
    int proc = atoi(argv[2]);
    
    std::cout << "Intializing " << n << " x " << n << " matrix of random doubles...\n";
    time_t ct1 = time(0);
    std::cout << "Clock time: " << ct1 << "\n";
    
    
    //double a[n][n];
    //i think im seg faulting form too much stack usage...also i now have to do 1d arrays since c sucks and wont let me do 2d dynamo arrays
    double *a = new double[n*n];
    double *a_orig = new double[n*n];
    
    omp_set_num_threads(proc);
    
    #pragma	omp	parallel
    {
        int	iii	=	omp_get_thread_num();
        printf("Hello	from	thread	%d\n",	iii);
    }
    
    for(int innitI = 0; innitI < n; innitI++)
    {
        for(int innitJ = 0; innitJ < n; innitJ++)
        {
            a[innitI*n+innitJ] = drand48();
            a_orig[innitI*n+innitJ] = a[innitI*n+innitJ];
        }
    }
    
    
    
    
    std::cout << "Done Initializing matrix a!\n";
    
    time_t ct2 = time(0);
    double ctd1 = difftime(ct1, ct2) * -1.0;
    std::cout << "Matrix a random initialization took: " << ctd1 << " second(s)\n";
    
    //std::cout << "Initializing a matrix took this much time: " << ctd1 << "\n";
    
    if(verbose)
    {
        std::cout.precision(3);
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << a[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
    }
    
    std::cout << "Begin running pseudo code...\n";
    double *pi = new double[n];
    double *l = new double[n*n];
    double *u = new double[n*n];
    
    //initializations
    for(int innitI = 0; innitI < n; innitI++)
    {
        for(int innitJ = 0; innitJ < n; innitJ++)
        {
            u[innitI*n+innitJ] = 0;
        }
    }
    for(int innitI = 0; innitI < n; innitI++)
    {
        for(int innitJ = 0; innitJ < n; innitJ++)
        {
                l[innitI*n+innitJ] = 0;

        }
    }
    if(verbose)
    {
        std::cout << "initialized u matrix: \n";
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << u[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
        
        std::cout << "initialized l matrix: \n";
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << l[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
    }
    
    
    //end of innitializations
    //std::cout << "WTF";
    //std::cout << "WTF" << std::flush;
    //begin pseudo
    int kPrime;
    
    for(int i = 0; i < n; i++)
    {
        pi[i] = i;
    }
    //std::cout << "seg fault";
    for(int k = 0; k < n; k++)
    {
        double max = 0;
        for(int i = k; i < n; i++)
        {
            //std::cout << "outter1 " << max << "\n";
            if(max < std::abs(a[i*n+k]))
           {
               
               max = std::abs(a[i*n+k]);
               kPrime = i;
               //std::cout << "inner1 " << max << "\n";
               //std::cout << "inner2 " << std::abs(a[i*n+k]) << "\n";
           }
        }
        if(max == 0)
        {
            std::cout << "ERROR (SINGULAR MATRIX)";
            //std::flush;
            //exit 0;
        }
        
        //swap π[k] and π[k']
        double swapMe = pi[k];
        pi[k] = pi[kPrime];
        pi[kPrime] = swapMe;
        
        int swapI;
        double swapMe2;
        double swapMe3;
        /*
        int chunk = ( n / proc );
        #pragma omp parallel shared(a,k,n,kPrime,chunk) private(swapI,swapMe2,swapMe3)
        {
            
        #pragma omp for schedule(dynamic,chunk)
         */
        for(swapI = 0; swapI < n; swapI++)
        {
            swapMe2 = a[k*n+swapI];
            a[k*n+swapI] = a[kPrime*n+swapI];
            a[kPrime*n+swapI] = swapMe2;
        }

        //#pragma omp for schedule(dynamic,chunk)
        for(swapI = 0; swapI < n; swapI++)
        {
            swapMe3 = l[k*n+swapI];
            l[k*n+swapI] = l[kPrime*n+swapI];
            l[kPrime*n+swapI] = swapMe3;
        }
            
        //}
        
        //continue past swaps
        u[k*n+k] = a[k*n+k];
        for(int i = k+1; i < n; i++)
        {
            l[i*n+k] = a[i*n+k] / u[k*n+k];
            u[k*n+i] = a[k*n+i];
        }
        
        int chunk = (( n - (k+1)) / proc) + 1; //adding plus one to ensure never zero and some procs will just have to figure out they ran out of work at the end
        //std::cout << "chunk size: " << chunk << "\n";
        int t_i;
        int t_j;
        #pragma omp parallel shared(a,k,n,l,u,chunk) private(t_i,t_j)
        //#pragma omp parallel shared(a,l,u,chunk) private(k,n,t_i,t_j)
        {
        #pragma omp for schedule(dynamic,chunk)
        for(t_i = k+1; t_i<n; t_i++)
        {
            for(t_j = k+1; t_j<n; t_j++)
            {
                a[t_i*n+t_j] = a[t_i*n+t_j] - (l[t_i*n+k] * u[k*n+t_j]);
            }
        }
            
        }
        
    }
    
    //put 1s on l's diag
    for(int innitI = 0; innitI < n; innitI++)
    {
        for(int innitJ = 0; innitJ < n; innitJ++)
        {
            if(innitI == innitJ)
            {
                l[innitI*n+innitJ] = 1;
            }else
            {
                //do nothing
            }
            
        }
    }
    
    
    time_t ct3 = time(0);
    double ctdAlgTotal = difftime(ct2, ct3) * -1.0;
    std::cout << "Algorithm in full took: " << ctdAlgTotal << " second(s)\n";
    
    
    double *pa = new double[n*n];
    double *lu = new double[n*n];
    
    //how do they look
    if(verbose)
    {
        std::cout << "Filled u matrix: \n";
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << u[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
        
        std::cout << "Filled l matrix: \n";
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << l[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
        
        std::cout << "filled pi vector/array: \n";
        for(int innitI = 0; innitI < n; innitI++)
        {
                std::cout << pi[innitI] << "\t";
        }
        std::cout << "\n";
    }
    
        //also fill the p matrix
        double *p = new double[n*n];
        for(int innitI = 0; innitI < n; innitI++)
        {
            int putOneHere = pi[innitI];
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                if(innitJ == putOneHere)
                {
                    p[innitI*n+innitJ] = 1;
                 //std::cout << "1\t";
                }else
                {
                    p[innitI*n+innitJ] = 0;
                  //std::cout << "0\t";
                }
            }
            //std::cout << "\n";
        }
    
    //how do they look
    if(verbose)
    {
        std::cout << "what pi looks like as a matrix: \n";
        //now actually show (just added calc)
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << p[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
    }
    
        std::cout << "\n";
    
    //how do they look
    if(verbose)
    {
        std::cout << "post operations a matrix(doesnt matter, only original matters really): \n";
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << a[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
    }
    
        //std::cout << "post operations a matrix using new func: \n";
        //printMatrix(a,n);
        
        
        
    
    
        if(showL21)
        {
        //calculate pa
        for(int leftRow = 0; leftRow < n; leftRow++)
        {
            for(int rightCol = 0; rightCol < n; rightCol++)
            {
                double total = 0.0;
                for(int eachElem = 0; eachElem < n; eachElem++)
                {
                    total = total + p[leftRow*n+eachElem] * a_orig[eachElem*n+rightCol];
                }
                pa[leftRow*n+rightCol] = total;
            }
        }
    
    //how do they look
    if(verbose)
    {
        std::cout << "Verification: \n";
        std::cout << "p * a: \n";
        //print pa
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << pa[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
    }
    
        
    
        //calculate lu
        for(int leftRow = 0; leftRow < n; leftRow++)
        {
            for(int rightCol = 0; rightCol < n; rightCol++)
            {
                double total = 0.0;
                for(int eachElem = 0; eachElem < n; eachElem++)
                {
                    total = total + l[leftRow*n+eachElem] * u[eachElem*n+rightCol];
                }
                lu[leftRow*n+rightCol] = total;
            }
        }
    
    //how do they look
    if(verbose)
    {
        std::cout << "l * u: \n";
        //print lu
        for(int innitI = 0; innitI < n; innitI++)
        {
            for(int innitJ = 0; innitJ < n; innitJ++)
            {
                std::cout << lu[innitI*n+innitJ] << "\t";
            }
            std::cout << "\n";
        }
    }
        }
    
    if(showL21)
    {
        double l21Norm = 0.0;
        for(int normCol = 0; normCol < n; normCol++)
        {
            double eucNorm = 0.0;
            for(int normRow = 0; normRow < n; normRow++)
            {
                eucNorm = eucNorm + pow((pa[normRow*n+normCol] - lu[normRow*n+normCol]),2.0);
            }
            
            eucNorm = pow(eucNorm,0.5);
            l21Norm = l21Norm + eucNorm;
        }
        
        std::cout.precision(7);
        std::cout << "The L21 Norm value of P*A - L*U is " << l21Norm << "\n";
        
    }
    time_t ct4 = time(0);
    double ctdVerif = difftime(ct3, ct4) * -1.0;
    std::cout << "ct3: " << ct3 << "\n";
    std::cout << "ct4: " << ct4 << "\n";
    std::cout << "Verification and l21 norm calc took: " << ctdVerif << " second(s)\n";
    
    
    //delete stuff
    delete [] a;
    delete [] a_orig;
    delete [] l;
    delete [] u;
    delete [] p;
    delete [] pa;
    delete [] lu;
    delete [] pi;
    
    //next?
    return 0;
}

