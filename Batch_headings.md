# Running jobs on the cluster

Please see [the Research Computing website](https://rc-docs.northeastern.edu/en/latest/runningjobs/index.html) for full information.

Before writing your batch script it's good to note the version numbers of any modules you'll ask it to load, since from within nano you won't be able to use tab-complete. Also, you do not need to be on a computing node to run the script, and you can specify any modules needed from within the script without loading them into your environment yourself.

## Headings for your batch script:

```
#!/bin/bash

#SBATCH –nodes=1
#SBATCH –time=4:00:00
#SBATCH –job-name=MyCPUJob
#SBATCH –partition=courses
#SBATCH –mail-type=ALL
#SBATCH –mail-users=username@northeastern.edu

# Load necessary modules, if needed

# Go to specific location, if needed

# Command
```
-J is the jobname; name it whatever you like. -N is the number of nodes you're requesting and -n is the number of threads. -o gives the form of the output file, which in this case will say "output_" followed by the job ID# assigned to the task. Likewise, -e gives the form of the output file.

After saving the script, run it with the command:

```sbatch script_name.sh```

Check on its status with the command:
```squeue -u user_name```
For example, I'd use ```squeue -u jhenzy```

When the job runs, output and error files will appear in the folder in which you ran the script. Check their contents for any error messages, or to see that everything is running according to plan.

If you need to cancel a job, use the command ```scancel <job_ID>```, supplying the specific job ID#.

