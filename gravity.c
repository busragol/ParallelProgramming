#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <time.h>

#define BILLION 1E9

int main(int argc, char *argv[])
{
    int processRank, processCount;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &processRank );
    MPI_Comm_size( MPI_COMM_WORLD, &processCount );

    int counter, counter2, starCount = 0, starsWanted = (int)(strtol(argv[argc - 1], NULL, 10));
    double R, G = 0.0006674, startTime, endTime, totalTime;

    struct timespec start, end;

    if(processRank == 0)
    {
        int xTemp, yTemp, zTemp, weightsTemp;

        FILE *inputFile;

        if( (inputFile = fopen(argv[argc - 2] , "r")) != NULL )
        {
            while(fscanf(inputFile, "%d\t%d\t%d\t%d", &xTemp, &yTemp, &zTemp, &weightsTemp) != -1)
            {
                starCount++;
            }

            rewind(inputFile);

            if(starsWanted > starCount)
            {
                printf("\"%s\" does not have more than %d stars.\n", argv[argc - 2], starCount);
                return(11);
            }
        }

        else
        {
            printf("The specified file \"%s\" does not exist.", argv[argc - 2]);
            return(1);
        }
    }

    MPI_Bcast(&starCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int xLocations[starsWanted], yLocations[starsWanted], zLocations[starsWanted], weights[starsWanted];
    double F[starsWanted];

    for(counter = 0 ; counter < starsWanted ; counter++)
    {
        xLocations[counter] = -1;
        yLocations[counter] = -1;
        zLocations[counter] = -1;
        weights[counter] = -1;
        F[counter] = 0.0;
    }

    if(processRank == 0)
    {
        FILE *inputFile;

        if( (inputFile = fopen(argv[argc - 2] , "r")) != NULL )
        {
            for(counter = 0 ; counter < starsWanted ; counter++)
            {
                fscanf(inputFile, "%d\t%d\t%d\t%d", &xLocations[counter], &yLocations[counter], &zLocations[counter], &weights[counter]);
            }
        }

        else
        {
            printf("The specified file \"%s\" does not exist.", argv[argc - 2]);
            return(1);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD );
	startTime = MPI_Wtime();

    MPI_Bcast(xLocations, starsWanted, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(yLocations, starsWanted, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(zLocations, starsWanted, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(weights, starsWanted, MPI_INT, 0, MPI_COMM_WORLD);

    /*for(counter = 0 ; counter < starsWanted ; counter++)
        {
            printf("core %d %d\t%d\t%d\t%d\n", rank, xLocations[counter], yLocations[counter], zLocations[counter], weights[counter]);
        }*/
    int loopStart, loopEnd;
    int remainder = starsWanted % processCount; //calculates how much data is left after an even distribution of work

    //adds one more variable to calculate to the list of variables that this group of threads need to calculate (distributes the leftover variables from even distribution)
    if(processRank < remainder) //these are threads that will get an extra variable
    {
        loopStart = processRank * ( (starsWanted / processCount) + 1 );

        loopEnd = (processRank + 1) * ( (starsWanted / processCount) + 1 );
    //printf("call sum, if 1, tp %d, ls %d, le %d, rm %d\n" , thread_part , loopStart , loopEnd , remainder);
    }

    else //the rest of the threads operate as usual with normal amount of data
    {
        loopStart = ( processRank * (starsWanted / processCount) ) + remainder;

        loopEnd = ( (processRank + 1) * (starsWanted / processCount) ) + remainder;
    //printf("call sum, if 2, tp %d, ls %d, le %d, rm %d\n" , thread_part , loopStart , loopEnd , remainder);
    }

    int localNumOfStars = loopEnd - loopStart;

    double localF[localNumOfStars]; //the force array of the calculated stars

    for(counter = 0 ; counter < localNumOfStars ; counter++)
    {
        localF[counter] = 0.0;
    }

    for(counter = loopStart ; counter < loopEnd ; counter++)
    {
        for(counter2 = 0 ; counter2 < starsWanted ; counter2++)
        {
            R = 0.0;

            if(counter == counter2)
            {
                continue;
            }

            else
            {
                R = sqrt(pow((xLocations[counter] - xLocations[counter2]), 2) + pow((yLocations[counter] - yLocations[counter2]), 2) + pow((zLocations[counter] - zLocations[counter2]), 2));

                localF[counter - loopStart] += ((G * weights[counter] * weights[counter2]) / pow(R, 2)); //make the index start from 0
                //printf("%lf\n",F[counter]);
            }
        }
    }

    /*for(counter = 0 ; counter < localNumOfStars ; counter++)
        {
            printf("%d %lf\n", processRank, localF[counter]);
        }*/

    MPI_Gather(localF, localNumOfStars, MPI_DOUBLE, F, localNumOfStars, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	endTime = MPI_Wtime();
	double timePassed = endTime - startTime;

	//printf("%d %lf\n", processRank, timePassed);

	MPI_Reduce(&timePassed, &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    /*if(processRank == 0)
    {
        for(counter = 0 ; counter < starsWanted ; counter++)
        {
            printf("%lf\n", F[counter]);
        }
    }*/

    if(processRank == 0)
    {
        printf("Average completion time of all processes: %lf seconds\n", totalTime / (double)processCount);

        FILE *outputFile;

        if( (outputFile = fopen("forces.txt" , "w")) != NULL )
        {
            for(counter = 0 ; counter < starsWanted ; counter++)
            {
                fprintf(outputFile, "%lf\n", F[counter]);
            }
        }

        else
        {
            printf("I/O error while saving forces.");
            return(2);
        }
    }

    MPI_Finalize();
}
