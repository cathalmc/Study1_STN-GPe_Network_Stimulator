import sys


newloc=sys.argv[1]
runtime=sys.argv[2]
cores=sys.argv[3]
net=sys.argv[4]


assert net in ["Small_world","Scale_free","Spatial","Regular","SBlock"]

fin = open("mRun.sh", "rt")
fout = open("Run.sh", "wt", newline='')

for line in fin:
    if line.strip()=="#SBATCH -n 10":
        new_core = "#SBATCH -n "+str(cores)
        fout.write(new_core+"\n")
    
    elif line.strip()=="#SBATCH -t TIME":
        new_time = "#SBATCH -t "+str(runtime)+":00:00"
        fout.write(new_time+"\n")
    elif line.strip()=="model_dir=\"${HOME}/storage/ClusterRuns\"":
        new_dir="model_dir=\"${HOME}/storage/"+str(newloc) +"\""
        fout.write(new_dir+"\n")
        
    elif line.strip()=="#SBATCH --job-name=SWS":
        new_name = "#SBATCH --job-name=" + net+"_"+ str(newloc)
        fout.write(new_name+"\n")
        
    elif line.strip()=="NETWORK=":
        newnetl = "NETWORK=" +"\"" +net  +"\""
        fout.write(newnetl+"\n")
        
    else:
        fout.write(line)

fin.close()
fout.close()

