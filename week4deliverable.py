# PRA2003 Week 4 Deliverable
# some corrections added from last week!!
# 1. reads each file, tallies every particle code per event, excluding empty events (header says 0 particles) from N.
# 2. combines the 10 files' results into one final average per code, using the sub-sampling method
#    (average of totals; uncertainty = std of the 10 per-file averages).
# 3. for every code / -code pair, tests whether the difference is significant, using the SAME
#    sub-sampling method on the difference itself (keeps the correlation between paired codes -
#    see "correlation_r" in the CSV - instead of combining each side's uncertainty independently,
#    which would be wrong if the two are correlated across files).
#
# Saves: subsample_results.csv, results.csv, significance.csv
# Full detail (every code, every pair) is in those CSVs - console output below is just a summary.
#
# NOTE: Make sure output-Set1.txt through output-Set10.txt are all in your directory.

import math        
import csv         # for writing the output tables to .csv files
import statistics  # for stdev()

# files & settings 
file_names = ["output-Set" + str(i) + ".txt" for i in range(1, 11)]  # output-Set1.txt ... output-Set10.txt
subsample_csv = "subsample_results.csv"   # per-file results (one row per code, per file)
results_csv = "results.csv"               # combined result: one row per code, averaged over all files
significance_csv = "significance.csv"     # pairwise asymmetry test results

threshold = 3   # (for later) a pair only counts as significant if its difference is more than
                # this many standard deviations (n_sigma) away from zero

progress_every = 500000   # print a progress line every this many events (per file),
                           

# known codes and their names 
# only these 12 codes are recorded
known_codes = [211, -211, 321, -321, 2212, -2212,
               3122, -3122, 3312, -3312, 3334, -3334]
known_names = ["Carbon monoxide", "Carbon-13 monoxide",
               "Nitric oxide", "Ionised NO",
               "Water", "Heavy water (D2O)",
               "Methane (CH4)", "Methyl ion (CH3-)",
               "Ethylene (C2H4)", "Ionised ethylene (C2H3-)",
               "Ozone (O3)", "Superoxide anion (O2-)"]
name_lookup = dict(zip(known_codes, known_names))   # code -> name, (e.g 211 -> carbon monoxide)
core_pairs = [211, 321, 2212, 3122, 3312, 3334]      # the 6 normal codes



# 1. read one file, tally every known code per (non-empty) event

#keep a running total per code:
   # total_count: the sum of counts, used to get the mean
   # Poisson: sqrt(N) / n_events

def analyse_file(filename):
    try:
        f = open(filename, "r") # protections!!
    except FileNotFoundError:
        print("Error: could not find the file", filename)
        return None   # it can be skipped

    # initialising
    total_count = {code: 0 for code in known_codes}
    n_events = 0   # only counts non-empty events 
    malformed_lines = 0   

    for line in f:   
        columns = line.split()

        if len(columns) == 2:   # a header should only have event id and particle count
            try:
                n_particles = int(columns[1]) # protection against non-integers
            except ValueError:
                print("Error, count is not an integer", line)
                malformed_lines += 1   # skip malformed header
                continue

            if n_particles > 0:   # empty event is not counted 
                n_events = n_events + 1 # tallying

                if n_events % progress_every == 0:   # checks that it is working 
                    print(f"  ... {filename}: {n_events:,} events read")

        elif len(columns) == 4:   # particle line usually has px, py, pz,  and code
            try:
                particle_id = int(columns[3])   # the code is the 4th column, python starts at 0 so 3 = 4
            except ValueError:
                malformed_lines += 1
                continue
            if particle_id in name_lookup:   # ignore unknown codes
                total_count[particle_id] += 1   

        else:   # skipped
            malformed_lines += 1

    f.close() # closing file

    if malformed_lines > 0:   
        print(f"Note: {malformed_lines} line(s) skipped in {filename}")

    if n_events == 0:   #if a file has only empty events 
        print("Error:", filename, "has no non-empty events so it is skipped")
        return None

    # running totals --> mean and an uncertainty for each code
    file_results = []
    for code in known_codes:
        mean = total_count[code] / n_events

        # Poisson: for a total count N, sigma = sqrt(N), so the uncertainty on the average is sqrt(N) / n_events   
        error = math.sqrt(total_count[code]) / n_events

        file_results.append({
            "code": code, "name": name_lookup[code],
            "n_events": n_events, "total_count": total_count[code],
            "average_per_event": float(f"{mean:.5g}"),   #  5 significant figures 
            "uncertainty": float(f"{error:.5g}"),        
            "subsample": filename
        })
    return file_results


# doing this on each file
all_results = []
for filename in file_names:
    file_results = analyse_file(filename)
    if file_results is not None:   # none means that file failed to read but is skipped
        all_results.extend(file_results)

if len(all_results) == 0:
    print("No files were read")
    exit()

# save the per-file results to a CSV 
# write turns dictionaries into CSV rows
with open(subsample_csv, "w", newline="") as out_f:
    writer = csv.DictWriter(out_f, fieldnames=["code", "name", "n_events", "total_count",
                                                "average_per_event", "uncertainty", "subsample"])
    writer.writeheader()
    writer.writerows(all_results)

files_read = len(set(row["subsample"] for row in all_results))   # counts only the files that were actually read
print(f"Part 1 done: {files_read} of {len(file_names)} files read -> {subsample_csv}") # to see whilst you are running the code!



# 2. combine the 10 files into one final average per code


# how many (non-empty) events were in each file for per-event average
events_per_file = {row["subsample"]: row["n_events"] for row in all_results} # builds a dictionary by looping over all_results
sub_sample_names = sorted(events_per_file.keys())   # sorted() puts them file names in a consistent order

# regroup the per-file rows for a dictionary of dictionaries

counts_by_code_file = {} # e.g  look up code X total count in output-SetX.txt 
names_by_code = {}
for row in all_results:
    counts_by_code_file.setdefault(row["code"], {})[row["subsample"]] = row["total_count"]
    names_by_code.setdefault(row["code"], row["name"])

# average per event in each file in order of subsamples
def per_file_average(code):
    file_counts = counts_by_code_file.get(code, {})
    return [file_counts.get(filename, 0) / events_per_file[filename] for filename in sub_sample_names]


# sub-sampling calculation 
final_results = []
for code in counts_by_code_file:
    file_counts = counts_by_code_file[code]

    # total count across files/divided by total events across files (pooled average)
    total_count = sum(file_counts.get(filename, 0) for filename in sub_sample_names)
    average_per_event = total_count / sum(events_per_file.values())

    # uncertainty is thestandard deviation of the 10 per-file averages
    # no assumption about correlation
    sub_sample_averages = per_file_average(code)
    uncertainty = statistics.stdev(sub_sample_averages) if len(sub_sample_averages) > 1 else 0.0

    final_results.append({
        "code": code, "name": names_by_code[code], "total_count": total_count,
        "average_per_event": round(average_per_event, 6), "uncertainty": round(uncertainty, 6)
    })

final_results.sort(key=lambda row: row["average_per_event"], reverse=True)   # most common first
final_by_code = {row["code"]: row for row in final_results}   # quick lookup by code for the last section of this code

with open(results_csv, "w", newline="") as out_f: # using the csv
    writer = csv.DictWriter(out_f, fieldnames=["code", "name", "total_count",
                                                "average_per_event", "uncertainty"])
    writer.writeheader() # writer again!
    writer.writerows(final_results)

print(f"Part 2 done: {sum(events_per_file.values()):,} events -> {results_csv}\n") # to keep track

# print the table of results to the console 
header_fmt = "{:>8} {:<30} {:>14} {:>18} {:>12}" # formatting
row_fmt = "{:>8} {:<30} {:>14} {:>18.6f} {:>12.6f}"
print(header_fmt.format("code", "name", "total_count", "average_per_event", "uncertainty"))
for row in final_results:
    print(row_fmt.format(row["code"], row["name"], row["total_count"],
                          row["average_per_event"], row["uncertainty"]))



# 3. asymmetry test

# Sample standard deviation 
# 0 if there's fewer than 2 
def safe_stdev(values):
    return statistics.stdev(values) if len(values) > 1 else 0.0

# Pearsons correlation coefficient between  two codes per file averages)
# +1 means they move together across files
# close to 0 means independent
# Why do we do this: tells us whether to combine their uncertainties separately, or whether we need the correlated method b
def pearson_r(x, y):
    mean_x, mean_y = sum(x) / len(x), sum(y) / len(y)
    cov = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y)) / (len(x) - 1)
    sd_x, sd_y = safe_stdev(x), safe_stdev(y)
    return cov / (sd_x * sd_y) if sd_x > 0 and sd_y > 0 else 0.0


# which codes have a pair
paired_codes = [code for code in final_by_code if code > 0 and -code in final_by_code]

# make sure none of the 6 molecules don't have a pair 
missing_pairs = [code for code in core_pairs if code not in paired_codes]
if missing_pairs:
    print("Error: missing partner for code(s):", missing_pairs)
    exit()

# put the 6 main pairs first, then any other paired codes sorted by how common they are
other_codes = sorted((c for c in paired_codes if c not in core_pairs),
                      key=lambda c: final_by_code[c]["total_count"], reverse=True)
pair_codes = core_pairs + other_codes

results = []
for code in pair_codes:
    a, b = final_by_code[code]["average_per_event"], final_by_code[-code]["average_per_event"]
    difference = a - b
    asymmetry = difference / (a + b)   # normalised difference (fraction of the total)

    # geach code's per-file averages, how the difference behaves across files directly 
  
    file_a, file_b = per_file_average(code), per_file_average(-code)
    per_file_diff = [fa - fb for fa, fb in zip(file_a, file_b)]
    per_file_asym = [(fa - fb) / (fa + fb) for fa, fb in zip(file_a, file_b) if (fa + fb) != 0]

    r = pearson_r(file_a, file_b)                # how correlated the pair is across files
    sigma = safe_stdev(per_file_diff)             # uncertainty on the difference
    sigma_asymmetry = safe_stdev(per_file_asym)   # uncertainty on the asymmetry

    # how many sigma the difference is away from zero for signficance
    n_sigma = abs(difference) / sigma if sigma > 0 else (float("inf") if difference != 0 else 0.0)

    # results is a list that started empty .append(...) adds one dictionary onto the end of that list
    results.append({
        "code": code, "name": final_by_code[code]["name"], "partner_name": final_by_code[-code]["name"],
        "avg_code": a, "avg_minus_code": b,
        "difference": round(difference, 4), "uncertainty": round(sigma, 3), "correlation_r": round(r, 3),
        "n_sigma": round(n_sigma, 2), "asymmetry_percent": round(100 * asymmetry, 3),
        "asymmetry_uncertainty_percent": round(100 * sigma_asymmetry, 3),
        "significant": n_sigma > threshold   # true if the difference is >3sigma
    })

with open(significance_csv, "w", newline="") as out_f: # now opening the significance.csv
    writer = csv.DictWriter(out_f, fieldnames=["code", "name", "partner_name",
                                                "avg_code", "avg_minus_code", "difference", "uncertainty",
                                                "correlation_r", "n_sigma", "asymmetry_percent",
                                                "asymmetry_uncertainty_percent", "significant"])
    writer.writeheader() # writer again again
    writer.writerows(results)

print(f"\nPart 3 done -> {significance_csv}\n") # checking if still running

# print one line per pair to not open CSV
for row in results:
    verdict = "significant" if row["significant"] else "not significant"
    print(f"{row['name']} (code {row['code']}): "
          f"difference = {row['difference']} +/- {row['uncertainty']}, "
          f"n_sigma = {row['n_sigma']}, r = {row['correlation_r']} -> {verdict}")

# summary in a sentence
n_sig = sum(row["significant"] for row in results)
sig_names = [row["name"] for row in results if row["significant"]]
print(f"\n{n_sig} of {len(results)} pairs significant (>{threshold} sigma):",
      ", ".join(sig_names) if sig_names else "none")
