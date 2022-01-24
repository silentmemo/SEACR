#### This script is based on the SEACR_1.3.sh. It is written in python form for my readability. This version is a lot more memory hungry, and thus should not be used for production. For Education use only. 

from json import load
from statistics import mode
import sys
def load_bedgraph_to_list(bdg_file_path):
    outList = []
    with open(bdg_file_path, mode="r")as inputFile:
        for line in inputFile:
            line = line.strip().split("\t")
            line[1] = int(line[1])
            line[2] = int(line[2])
            line[3] = float(line[3])
            outList.append(line)
    return outList
def make_signal_block(bedgraph_list,filename):
    ### this function make signal block, it joins consecutive record, record the max with in the block. count the number of consecutive block, calculate the cumulative signal (AUC), then print a bed file. Note that loading large list into memory is not wise. 
    print("Generating signal block file from: ", filename)
    done_flag = False
    record_visited = 0
    number_of_record = len(bedgraph_list)
    outList = []
    while not done_flag:
        
        ### set current record
        if record_visited >= number_of_record:
            done_flag = True
            break
        else:
            current_record =  bedgraph_list[record_visited]
            current_seqname = current_record[0]
            current_start = current_record[1]
            current_end = current_record[2]
            current_signal = current_record[3]

            current_AUC = current_signal * (current_end-current_start)
            current_num = 1 ### only visited one record 
            current_max = current_signal

            max_start = current_start
            max_end = current_end
            max_coordinate = current_seqname + ":" + str(max_start) + "-" + str(max_end)
            record_visited = record_visited +1

            for fwd_record in bedgraph_list[record_visited:]:
                fwd_seqname = fwd_record[0]
                fwd_start = fwd_record[1]
                fwd_end = fwd_record[2]
                fwd_signal = fwd_record[3]
                
                if fwd_start == current_end and fwd_seqname == current_seqname:
                    ### is consecutive
                    current_AUC = current_AUC + fwd_signal *(fwd_end-fwd_start) ### becasue it is cumulative, it update the current_AUC variable
                    
                    current_num = current_num +1
                    if fwd_signal > current_max:
                        max_start = fwd_start
                        max_end = fwd_end
                        max_coordinate = fwd_seqname + ":" + str(max_start) + "-" + str(max_end)
                        current_max = fwd_signal
                    elif fwd_signal == current_max:
                        max_end = fwd_end ### only update the max_end but exclude the max_start 
                        max_coordinate = fwd_seqname + ":" + str(max_start) + "-" + str(max_end)
                    current_end = fwd_end 
                    current_seqname = fwd_seqname
                    record_visited = record_visited +1 
                else:
                    break ### no need to keep looking if it is not consecutive, therefore bdg assume coordinate sorted
            ### output the signal blocks for this consecutive signal block 
            outList.append([current_seqname,str(current_start),str(current_end),str(current_AUC),str(current_max),str(max_coordinate),str(current_num)])
    return outList
def make_AUC_list(signal_block_list,filename):
    print("Generating .auc file from: ", filename)
    outList = []
    for record in signal_block_list:
        outList.append([float(record[3]), int(record[6])])
    return outList
def calculate_threshold_with_normalized_control():
    return None ### TODO 
def main():
    print("Running python version of SEACR_1.3.sh")
    target_bdg_file = sys.argv[1]
    target_list = load_bedgraph_to_list(target_bdg_file)
    target_signal_block_list = make_signal_block(target_list, "target")

    control_bdg_file = sys.argv[2]
    control_list = load_bedgraph_to_list(control_bdg_file)
    control_signal_block_list = make_signal_block(control_list, "control")

    with open("target.auc.bed",mode = "w")as outFile:
        for line in target_signal_block_list:
            outLine = "\t".join(line) + "\n"
            outFile.write(outLine)
    with open("control.auc.bed",mode = "w")as outFile:
        for line in control_signal_block_list:
            outLine = "\t".join(line) + "\n"
            outFile.write(outLine)
    target_auc_list = make_AUC_list(target_signal_block_list,"target")
    control_auc_list = make_AUC_list(control_signal_block_list,"control")

    
main()