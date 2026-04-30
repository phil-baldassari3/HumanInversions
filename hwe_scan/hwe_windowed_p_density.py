from collections import deque


### HELPER FUNCTIONS ###

def _joinany(sep, mylist):
    """Function to join a list by a sep regardless of element types"""

    strs = [str(x) for x in mylist]

    return sep.join(strs)



def _pop_bp_step(d, step):
    """
    Helper function used to pop a number of elements from the left
    of a deque equal to the step size in. This is done in place.
    """

    #count how many elements to pop
    count_ele_to_pop = 0
    for idx in range(len(d)):
        if d[idx].position <= step:
            count_ele_to_pop += 1
        else:
            break
    
    #pop from left
    for i in range(count_ele_to_pop):
        d.popleft()


def _compute_sigp_density(deque_window, thresh):
    """
    Helper function that computes the proportion of p-values in a deque
    that are below a specified threshold
    """

    #count significant p values
    ps_below = [x.pval for x in deque_window if x.pval < thresh]
    sigPs = len(ps_below)

    #count total p values
    totPs = len(deque_window)

    #significant p-value density
    den = sigPs/totPs

    return den


### CLASS ###

class HWELine():
    def __init__(self, hwe_line):
        """
        This class is for parsing a line from the plink --hardy output.

        :param hwe_line: record line from vcf
        :type hwe_line: str

        Attributes:
            chromosome (int): record chromosome
            position (int): record position
            pval (float): record HWE Chi-square p-value
        """

        self.chromosome = int(hwe_line.strip().split()[0])
        self.position = int(hwe_line.strip().split()[1])
        self.pval = float(hwe_line.strip().split()[-1])




### SLIDING WINDOW SCAN FUNCTION ###

def windowed_p_density_scan(plink_hwe, window_size, window_step, significance_threshold):
    """
    Function takes the plink --hardy output and computes the proportion of sites with a p-value
    below a given threshold (density of significant p-values per window). Note the sites much be
    sorted.
    
    :param vcfgz: gziped vcf file name or path to file.
    :param window_size: number of bps defining the window size.
    :param window_step: number of bps to slide window, recommended to be 10% of window_size.
    :param significance_threshold: value below which p-values are considered "significant"

    :type vcfgz: str
    :type window_size: int
    :type window_step: int
    :type significance_threshold: float
    """

    #setting up deque for sliding window
    window_deque = deque()

    #open files
    with open(plink_hwe, "r") as vcf, open(f"hwe_{significance_threshold}pvalue_windowed_density_{window_size}bpwin_{window_step}step.bedgraph", "w") as outbedgraph:

        #window book keeping
        chrom = 1
        win_start = 1
        win_end = window_size
        win_step = window_step

        #loop through lines
        for line in vcf:
            #skip header lines
            if line.startswith(" CHR"):
                continue

            ###sliding window scan###

            #holding on  current record
            current_record = HWELine(line)

            #Check if we are still on the same chromsome
            if current_record.chromosome != chrom:

                #compute p density
                sigPdensity = _compute_sigp_density(window_deque, significance_threshold)
                #write to file
                outbedgraph.write(_joinany("\t", [f"chr{chrom}", win_start, win_end, sigPdensity]) + "\n")
                #empty deque
                window_deque = deque()
                #reset window bookeeping
                print(f"Done running through Chromosome {chrom}")
                chrom = current_record.chromosome
                win_start = 1
                win_end = window_size
                win_step = window_step

            else:

                if current_record.position < win_start:
                    exit("Record position before window start position")
                
                if current_record.position <= win_end:
                    window_deque.append(current_record)
                    continue

            #increment window until position is reached
            for i in range(win_end, current_record.position, window_step):

                if len(window_deque) == 0:
                    #write to file
                    outbedgraph.write(_joinany("\t", [f"chr{chrom}", win_start, win_end, 0]) + "\n")
                    
                else:
                    #compute p density
                    sigPdensity = _compute_sigp_density(window_deque, significance_threshold)
                    #sliding to next step
                    _pop_bp_step(window_deque, win_step)
                    #write to file
                    outbedgraph.write(_joinany("\t", [f"chr{chrom}", win_start, win_end, sigPdensity]) + "\n")

                #increment window
                win_end = i + window_step
                win_start = win_end - window_size + 1
                win_step = win_start + window_step
                
            #append current record to deque
            window_deque.append(current_record)

        #clean up last window
        if len(window_deque) > 0:
            #compute p density
            sigPdensity = _compute_sigp_density(window_deque, significance_threshold)
            #write to file
            outbedgraph.write(_joinany("\t", [f"chr{chrom}", win_start, win_end, sigPdensity]) + "\n")





def main():
    windowed_p_density_scan("hwe_persite.hwe", 100000, 10000, 0.01)


if __name__ == "__main__":
    main()
    print("DONE!")

