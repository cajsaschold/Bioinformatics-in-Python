import sys
import re
import matplotlib.pyplot as plot


def allele_freq(line, dict_pos_to_freq):
    columns = line.split("\t")
    info_column = columns[7]
    pos = columns[1]
    match = re.search('AF=([^;]+)', info_column)

    if match:
        freqValue = float(match.group(1)) #(1) to not include "AF=", converts from string to float
                
        dict_pos_to_freq[pos]=freqValue

    return dict_pos_to_freq


def freq_histo(dict_pos_to_freq):
    freqs = dict_pos_to_freq.values()

    plot.hist(freqs, bins=100, range=(0,1), log=True)
    plot.xlabel("Alternative Allele Frequency")
    plot.ylabel("Number of Variants")
    plot.title("Histogram for the Alternative Allele Frequency")
    plot.savefig('Freq hisogram')
    plot.clf() 

def x_per_sample(line, sample_names, dict_sample_to_var_counts, dict_sample_to_sing_counts, dict_sample_to_com_counts):
    columns = line.split("\t")
    info_column = columns[7]
    genotypes = columns[9:]
    match = re.search('AF=([^;]+)', info_column)
    if match:
        freqValue = float(match.group(1))

    for idx, genotype in enumerate(genotypes):
        if genotype in ["1|0", "0|1", "1|1"]:
            sample_name = sample_names[idx]
            dict_sample_to_var_counts[sample_name] += 1
            if 'AC=1;' in info_column:
                sample_name = sample_names[idx]
                dict_sample_to_sing_counts[sample_name] += 1
            if freqValue > 0.05:
                sample_name = sample_names[idx]
                dict_sample_to_com_counts[sample_name] += 1

    return dict_sample_to_var_counts, dict_sample_to_sing_counts, dict_sample_to_com_counts


def x_per_sample_plot(dict_sample_to_x, dict_sample_to_country, dict_country_to_continent, title):

    #creates two dicts where the key is the keys from the population-dict and each value is initialized to 0
    dict_country_to_x_total = {country: 0 for country in dict_country_to_continent.keys()}
    dict_country_to_sample_count= {country: 0 for country in dict_country_to_continent.keys()} #count of samples for each country.

    for sample, x_count in dict_sample_to_x.items():
        country = dict_sample_to_country[sample]
        dict_country_to_x_total[country] += x_count
        dict_country_to_sample_count[country] += 1

    #new dict by looping through all populations and setting the population as key and the avrage as value 
    dict_country_to_average= {country: dict_country_to_x_total[country] / dict_country_to_sample_count[country] for country in dict_country_to_x_total}

    sorted_dict = dict(sorted(dict_country_to_average.items(), key=lambda item: item[1]))
    plot.bar(sorted_dict.keys(), sorted_dict.values())
    plot.xticks(rotation='vertical')
    plot.ylabel("Average Number of " + title + " per Sample")
    plot.title("Average Number of " + title + " per Sample for Each Population")

    plot.savefig(title)
    plot.clf()


def count(line, ex_count, my_gene_count):
    columns = line.split("\t")
    info_column = columns[7]

    if 'Func.refGene=nonsynonymous_SNV' in info_column:
        ex_count += 1 
        if 'Gene.refGene=PTPN22' in info_column:
            my_gene_count += 1 

    return ex_count, my_gene_count


def read_vcf(filename):
    dict_pos_to_freq = {}
    ex_count = 0
    my_gene_count = 0

    with open(filename, 'r') as vcf:
        for line in vcf:
            line = line.strip()
        
            if line.startswith('#CHROM'):
                sample_names = line.split("\t")[9:]
                dict_sample_to_var_count = {sample: 0 for sample in sample_names} 
                dict_sample_to_sing_count = {sample: 0 for sample in sample_names}
                dict_sample_to_com_count = {sample: 0 for sample in sample_names}

            elif not line.startswith('#'):
                allele_freq(line, dict_pos_to_freq)

                ex_count, my_gene_count = count(line, ex_count, my_gene_count)
                
                x_per_sample(line, sample_names, dict_sample_to_var_count, dict_sample_to_sing_count, dict_sample_to_com_count)

    return dict_pos_to_freq, ex_count, my_gene_count, dict_sample_to_var_count, dict_sample_to_sing_count, dict_sample_to_com_count


def read1000(filename):
    sample = {}
    population = {}
    with open(filename, 'r') as sampleFile:
        for line in sampleFile:
            line = line.strip()
            if not line.startswith('Sample'):
                info = line.split('\t')
                sample[info[0]] = info[3]
                if info[3] not in population:
                    population[info[3]] = info[5]
    return sample, population


def main():
    CDS = 2424

    dict_pos_to_freq, ex_count, my_gene_count, dict_sample_to_var, dict_sample_to_sing, dict_sample_to_com = read_vcf(sys.argv[1])

    dict_sample_to_country, dict_country_to_continent = read1000(sys.argv[2])

    print("Total nr of samples:", len(dict_sample_to_var))
    print("Number of variants:", len(dict_pos_to_freq))
    print("Nonsynonymous variants:", ex_count)
    print("Nonsynonymous variantsin my gene:", my_gene_count)
    print("Nonsynonymous variants in my gene normalized by gene lenght:", my_gene_count/CDS)

    freq_histo(dict_pos_to_freq)

    x_per_sample_plot(dict_sample_to_var, dict_sample_to_country, dict_country_to_continent, 'Variants')
    x_per_sample_plot(dict_sample_to_sing, dict_sample_to_country, dict_country_to_continent, 'Singletons')
    x_per_sample_plot(dict_sample_to_com, dict_sample_to_country, dict_country_to_continent, 'Common Variants')

    

if __name__ == '__main__':
  main()
