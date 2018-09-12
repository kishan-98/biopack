# Usage: python script_biopac.py <evaluating python file> <source directory> <destination directory>
import multiprocessing as mp
import subprocess as sp
import random, string, math, sys, os, re, time, heapq

def analyse_text(eval_file, source_filelist, destination_filelist, idx):
    log = []
    i = 0
    for s_file, d_file in zip(source_filelist, destination_filelist):
        command = "python " + re.escape(eval_file) + " " + re.escape(s_file) + " " + re.escape(os.path.dirname(d_file))
        # command = "python gibber"
        # print("command:", command)
        try:
            exec_value = os.system(command)
            if exec_value:
                print(sp.check_output(command))
            else:
                # output.put((idx + i, s_file + ": Done"))
                log.append((idx + i, s_file + ": Done"))
        except:
            # output.put((idx + i, s_file + ": Error"))
            log.append((idx + i, s_file + ": Error"))
        i = i + 1
        print(s_file)
    return log

def k_partition(combined_filelist, k):
    source_filelist, destination_filelist, heap_filesize = [], [], []
    for i in range(k):
        heapq.heappush(heap_filesize, (0, i))
        source_filelist.append([])
        destination_filelist.append([])
    combined_filelist.sort(key=lambda detail: detail[2], reverse=True)
    for s_file, d_file, size_file in combined_filelist:
        heap_filesize_top = heapq.heappop(heap_filesize)
        source_filelist[heap_filesize_top[1]].append(s_file)
        destination_filelist[heap_filesize_top[1]].append(d_file)
        heapq.heappush(heap_filesize, (heap_filesize_top[0] + size_file, heap_filesize_top[1]))
    # print("heap_filesize:", heap_filesize)
    return source_filelist, destination_filelist

if __name__ == '__main__':
    start_time = time.time()

    output = mp.Queue()

    program_name = sys.argv[0]
    eval_file = ""
    source_directory = ""
    destination_directory = ""
    suffix = "_analysis"
    separater = "."
    extension = "txt"

    if len(sys.argv) < 2:
        print("No evaluating python file provided")
        sys.exit()
    elif len(sys.argv) < 3:
        print("No source directory provided")
        sys.exit()
    elif len(sys.argv) < 4:
        print("No destination directory provided; setting source directory as destination directory")
        eval_file = sys.argv[1]
        source_directory = destination_directory = sys.argv[2]
    else:
        eval_file = sys.argv[1]
        source_directory = sys.argv[2]
        destination_directory = sys.argv[3]

    source_filelist = []
    destination_filelist = []
    combined_filelist = []

    # Listing all the files
    for dir in os.listdir(source_directory):
        if not dir.startswith('.') and not os.path.isfile(os.path.join(source_directory, dir)):
            # print("dir:", dir)
            directory = os.path.join(destination_directory, dir + suffix)
            try:
                os.stat(directory)
            except:
                os.mkdir(directory)
            for file in os.listdir(os.path.join(source_directory, dir)):
                # print("file:", file)
                if file.startswith("crt_") and file.endswith(".txt"):
                    # print(os.path.abspath(os.path.join(os.path.join(source_directory, dir), file)))
                    source_file = os.path.abspath(os.path.join(os.path.join(source_directory, dir), file))
                    # source_file = re.escape(source_file)
                    # source_filelist.append(source_file)
                    destination_file = os.path.abspath(os.path.join(os.path.join(destination_directory, dir + suffix), os.path.splitext(os.path.basename(file))[0] + suffix + separater + extension))
                    # destination_file = re.escape(destination_file)
                    # destination_filelist.append(destination_file)
                    combined_filelist.append((source_file, destination_file, os.path.getsize(source_file)))

    # Sort all files based on the file size
    combined_filelist.sort(key=lambda filename: filename[2], reverse=True)

    # for s_file, d_file, filesize in combined_filelist:
        # print(s_file, d_file, filesize)
        # print(s_file, filesize)
        # print(filesize)
    # sys.exit()

    n = min(mp.cpu_count(), len(combined_filelist)) # # of partition to be made of filelist for optimized results

    # Rearranging files accordingly

    source_filelist, destination_filelist = k_partition(combined_filelist, n)

    # for i in range(n):
    #     source_filelist.append(list(map(lambda x: x[0], combined_filelist[i::n])))
    #     destination_filelist.append(list(map(lambda x: x[1], combined_filelist[i::n])))
    #     print(i, sum(list(map(lambda x: x[2], combined_filelist[i::n]))))

    # print(source_filelist)
    # print(destination_filelist)
    # print('\n\n'.join('\n'.join(s_list) for s_list in source_filelist))
    # print('\n\n'.join('\n'.join(d_list) for d_list in destination_filelist))
    # sys.exit()

    # processes = [mp.Process(target = analyse_text, args = (eval_file, source_filelist[x], destination_filelist[x], x, output)) for x in range(n)]
    #
    # for p in processes:
    #     p.start()
    #
    # for p in processes:
    #     p.join()
    #
    # results = [output.get() for p in processes]
    # results.sort()
    # results = [r[1] for r in results]
    # print('\n'.join(results))

    pool = mp.Pool(processes=n)
    output = [pool.apply_async(analyse_text, args = (eval_file, source_filelist[x], destination_filelist[x], x*math.ceil(len(source_filelist)/n),)) for x in range(n)]

    # output = [mp.Process(target = analyse_text, args = (eval_file, source_filelist[x*math.ceil(len(source_filelist)/n):(x+1)*math.ceil(len(source_filelist)/n)], destination_filelist[x*math.ceil(len(destination_filelist)/n):(x+1)*math.ceil(len(destination_filelist)/n)], x*math.ceil(len(source_filelist)/n),)) for x in range(n)]
    # for p in output:
    #     p.start()
    # for p in output:
    #     p.join()

    results = []
    for p in output:
        results = results + p.get()
    results.sort()
    results = [r[1] for r in results]
    print('\n'.join(results))
    end_time = time.time()

    print("Time taken:\t", math.floor((end_time - start_time)/60), "min", math.floor((end_time - start_time)%60), "sec")
