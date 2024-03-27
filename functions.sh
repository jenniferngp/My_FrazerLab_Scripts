print_exec() {
        local cmd=$1; shift
        local log=$1
        #local hold_jobs=$@ # save proceeding arguments
        #local njobs=0 # set njobs 
        date >& 2
        date >> $log
        echo $cmd >& 2
        echo $cmd >> $log
        eval $cmd
}
