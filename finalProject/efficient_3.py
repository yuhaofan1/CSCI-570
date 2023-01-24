import time
import psutil
import sys

delta_e = 30

alpha = {
    'A': {
        'A': 0,
        'C': 110,
        'G': 48,
        'T': 94
    },
    'C': {
        'A': 110,
        'C': 0,
        'G': 118,
        'T': 48
    },
    'G': {
        'A': 48,
        'C': 118,
        'G': 0,
        'T': 110
    },
    'T': {
        'A': 94,
        'C': 48,
        'G': 110,
        'T': 0
    }
}


def penalty(str1, str2):
    """calculate and print the penalty of two given string"""
    test_res = 0
    for i in range(len(str1)):
        if str1[i] == '_' or str2[i] == '_':
            test_res += delta_e
        else:
            test_res += alpha[str1[i]][str2[i]]
    print(test_res)


def process_memory():
    process = psutil.Process()
    mem_info = process.memory_info()
    mem = mem_info.rss/(1024)
    return mem

def string_generator(base_str: str, indices: list) -> str:
    res_str = base_str
    for index in indices:
        res_str = res_str[:index+1]+res_str+res_str[index+1:]
    return res_str


def test_string_generator():
    res_str = string_generator("ACTG", [3, 6, 1])
    true_res = "ACACTGACTACTGACTGGTGACTACTGACTGG"
    if res_str == true_res:
        print("Test success")
    else:
        print("Test wrong")


def read_param(filepath: str) -> dict:
    """read two strings' parameter from given file"""
    try:
        with open(filepath, "r") as fp:
            params = []
            param = {}
            base_str = ""
            for line in fp.readlines():
                line = line.strip()
                if line.isalpha():
                    if len(param) != 0:
                        params.append(param)
                        param={}
                    param[line] = []
                    base_str = line
                elif line.isdigit():
                    param[base_str].append(int(line))
            params.append(param)
            return params
    except FileNotFoundError:
        print("No that input file "+filepath)
        exit(0)


def test_read_param():
    filepath = "./UploadedProject/SampleTestCases/input1.txt"
    params = read_param(filepath)
    for k, v in params.items():
        print(string_generator(k, v))


def string_input() -> list:
    """generate a list from params contain two string as algorithm input """
    strings = []
    try:
        in_path = sys.argv[1]
        # in_path = "input4.txt"
    except IndexError:
        print("No input of filename")
        exit(0)
    a = {}
    params = read_param(in_path)
    for param in params:
        for k, v in param.items():
            strings.append(string_generator(k, v))
    return strings


def dp_last_column(seq1, seq2):
    """calculate the last column of seq alignment dp value"""
    l1 = len(seq1)
    l2 = len(seq2)
    dp = [0 for i in range(l2+1)]
    for i in range(1, l2+1):
        dp[i] = i*delta_e
    for i in range(1, l1+1):
        dpl = [dp[i] for i in range(len(dp))]
        dp[0] = i*delta_e
        for j in range(1, l2+1):
            dp[j] = min(dpl[j]+delta_e, dp[j-1]+delta_e,
                        dpl[j-1]+alpha[seq1[i-1]][seq2[j-1]])
    return dp


def find_opt_cut(str1, str2, xl, xr, yl, yr):
    """return the optimal slice index for str1 and str2"""
    mid = int((xl+xr)/2)
    dp_left = dp_last_column(str1[xl:mid+1], str2[yl:yr+1])
    dp_right = dp_last_column(str1[mid+1:xr+1][::-1], str2[yl:yr+1][::-1])
    min_cost = dp_left[0]+dp_right[-1]
    cur_cost = min_cost
    fh_letter_num = 0
    for i in range(1, len(dp_left)):
        cur_cost = dp_left[i]+dp_right[i*-1-1]
        if cur_cost < min_cost:
            min_cost = cur_cost
            fh_letter_num = i
    return mid, yl+fh_letter_num-1

def basic_seq_alignment(str1, str2):
    # initialization
    m = len(str1)
    n = len(str2)
    dp = [[0 for i in range(n+1)] for j in range(m+1)]
    # initialize dp table
    for i in range(m+1):
        dp[i][0] = delta_e*i
    for i in range(n+1):
        dp[0][i] = delta_e*i
    for i in range(1, m+1):
        for j in range(1, n+1):
            # ith letter's index is i-1
            dp[i][j] = min(dp[i-1][j-1]+alpha[str1[i-1]][str2[j-1]],
                           dp[i-1][j]+delta_e, dp[i][j-1]+delta_e)
    # Reconstructing alignment string
    str1_ali = ""
    str2_ali = ""
    i = m
    j = n
    while i > 0 and j > 0:
        if dp[i][j] == dp[i-1][j-1]+alpha[str1[i-1]][str2[j-1]]:
            str1_ali = str1[i-1]+str1_ali
            str2_ali = str2[j-1]+str2_ali
            i = i-1
            j = j-1
        elif dp[i][j] == dp[i-1][j]+delta_e:
            str1_ali = str1[i-1]+str1_ali
            str2_ali = '_'+str2_ali
            i = i-1
        else:
            str1_ali = '_'+str1_ali
            str2_ali = str2[j-1]+str2_ali
            j = j-1
    while i > 0:
        str1_ali = str1[i-1]+str1_ali
        str2_ali = '_'+str2_ali
        i = i-1
    while j > 0:
        str1_ali = '_'+str1_ali
        str2_ali = str2[j-1]+str2_ali
        j = j-1
    return str1_ali, str2_ali, dp[m][n]


def mesa(str1, str2, xl, xr, yl, yr):
    """
    divide and conquer method
    memory efficient sequence alignment
    """
    # Todo: base case
    if xr-xl < 2 or yr-yl < 2:
        return basic_seq_alignment(str1[xl:xr+1], str2[yl:yr+1])
    # normal case
    str1_cut_idx, str2_cut_idx = find_opt_cut(str1, str2, xl, xr, yl, yr)
    str1_fh_ali, str2_fh_ali, fh_cost = mesa(
        str1, str2, xl, str1_cut_idx, yl, str2_cut_idx)
    str1_sh_ali, str2_sh_ali, sh_cost = mesa(
        str1, str2, str1_cut_idx+1, xr, str2_cut_idx+1, yr)
    return str1_fh_ali+str1_sh_ali, str2_fh_ali+str2_sh_ali, fh_cost+sh_cost


def main():
    strings = string_input()
    str1 = strings[0]
    str2 = strings[1]
    l1 = len(str1)
    l2 = len(str2)
    # timing start
    start_time = time.time()
    # record memery usage start
    mem_before = process_memory()
    # Todo: sequence alignment algorithm
    str1_ali, str2_ali, cost = mesa(str1, str2, 0, l1-1, 0, l2-1)
    # record memory usage end
    mem_after = process_memory()
    # timing end
    end_time = time.time()
    # write into output file
    with open(sys.argv[2], "w") as fp:
        fp.write(str(cost)+'\n')
        fp.write(str1_ali+'\n')
        fp.write(str2_ali+'\n')
        fp.write(str((end_time-start_time)*1000)+'\n')
        fp.write(str(mem_after-mem_before)+'\n')


if __name__ == "__main__":
    main()
