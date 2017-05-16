# pp-matrix-multiplication

Use C and MPI to compute matrix multiplication A**k on a cluster of nodes, where A is an nXn matrix and k is a constantMatrix multiplication is one of the most used algorithm in scientific problems. So it’s important that we
devise an efficient parallel algorithm to perform this operation. So in this project we experimented with
two approaches and saw which is better of the both.

First the traditional 1-D block algorithm which has the worst time complexity then the cannon’s algorithm which is fast when compared to the 1D block algorithm.


**RUN:**

**_.pbs_** files are the job files that are submitted to the cluster where we specify the number of nodes, name of the executable file, and all other configurations for the job to run.

**_.pl_** files are the script files that automate the whole process by executing the c program and submitting the .pbs job to the cluster.

**_.err_** files give any error after executing the job.

**_.out_** files are the output files after the job being executed.

Run the .pl file to start the program execution on the server.
