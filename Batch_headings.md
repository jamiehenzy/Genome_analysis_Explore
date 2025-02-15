# Running jobs on the cluster

Please see [the Research Computing website](https://rc-docs.northeastern.edu/en/latest/runningjobs/index.html) for full information.

Before writing your batch script it's good to note the version numbers of any modules you'll ask it to load, since from within nano you won't be able to use tab-complete.

## Headings for your batch script:

```
#!/bin/bash
#SBATCH -J myjob                            # Job name
#SBATCH -N 2                                # Number of nodes
#SBATCH -n 16                               # Number of tasks
#SBATCH -o output_%j.txt                    # Standard output file
#SBATCH -e error_%j.txt                     # Standard error file

# Load necessary modules, if needed

# Go to specific location, if needed

# Command
```
-J is the jobname; name it whatever you like. -N is the number of nodes you're requesting and -n is the number of threads. -o gives the form of the output file, which in this case will say "output_" followed by the job ID# assigned to the task. Likewise, -e gives the form of the output file.

After saving the script, run it with the command:

``` sbatch script_name.sh ```

