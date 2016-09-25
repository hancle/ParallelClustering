# Parallel Clustering Algorithm 

Version 3 is prefered. 

Algorithm:

1. Broadcast array x_axis and y_axis to each processor. (provided by original codes).
2. Use max and min function to find the boundary of matrix.
3. Each core sorts the corresponding segment of array x_axis and y_axis, and use
MPI_Gather to collect them to p0. Find the midpoint.
4. Find the four corner coordinates of partitions based on the midpoint and boundary.
5. Every time when finish to get 2^i partitions, we broadcast the sorted array x_axis ad y_axis
to each processor in order to determine the direction (along x axis or y axis) of partition in
next iteration.
6. After we have the aimed number of quadrants, we record the boundary coordinates for each
quadrants and calculate the cost of each partition.


Scalability:

Our algorithm has proper scalability to fit different quadrants with different numbers of core. As shown in running time table, generally, implement with 32 processor gives highest performance.


Test Command:
example : mpiexec -n 4 /home/wyv390/358/bisection 64


Some explanations on codes:
In this project, we have 3 different versions based on the total cost calculation: (1) multi-core used to calculate the cost based on euclidean distance
(2) single-core used to calculate the cost based on euclidean distance
(3) multi-core used to calculate the cost based on the manhattan distance

The problem on the first version(multi-core, euclidean) is that, during the multiplication operations, it’s pretty easy to get overflow. That’s because the range for numbers in array X_axis and Y_axis (unsigned int data type) have been 2^10, if we implement operations like multiplication, the range will approach a very hight value such as 2^20. Furthermore, we also need to add those euclidean distance together to find the total cost. Since there are more than 520,000 numbers in original matrix, the total cost will be huge, which easy to cause the overflow. If we increase the range of data type(e.g. long int, long double), the running time will be super long. Due to this problem (without using long data type), by using algorithm version 1, we received different total costs with same number of quadrants and different number of cores.

The problem of version 2 is, it may also not be able to get correct values. At the same time, calculate the total cost with only one core will spend a lot of time.

Finally in version 3(preferred version), we would prefer to make a small change on the requirement: using manhattan distance instead of euclidean distance. Functionally, manhattan distance in this project can perform the same thing with euclidean— expressing the total cost. Same with euclidean distance, by using manhattan distance, we can observe the difference on cost between different quadrants and also compare the total costs on implementing different number of quadrants. Another benefit from manhattan distance is shorter running time compared with euclidean distance.