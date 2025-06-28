import streamlit as st
import matplotlib.pyplot as plt

st.title("DNA Sequence Analysis")


def read_fatsa(file):
    sequence = ""
    for line in file:
        line = line.decode("utf-8")
        if line.startswith(">"):
            continue
        sequence += line.strip()
    return sequence


userfile = st.file_uploader("Upload a FASTA file for analysis",
                            type=["txt", "fasta"])

if userfile:
    try:
        dna = read_fatsa(userfile)
        st.success("Successfully loaded DNA sequence.")
        st.write("DNA sequence:", dna)

        # Nucleotide count
        st.subheader("Nucleotide Analysis")
        d = {}
        for i in range(len(dna)):
            if dna[i] in d:
                d[dna[i]] += 1
            else:
                d[dna[i]] = 1
        st.write("Nucleotide count:", d)

        # Purine and Pyrimidine count
        totalpurinecount = 0
        totalpyrimidinecount = 0
        for item in d:
            if item == 'A' or item == 'G':
                totalpurinecount += d[item]
            elif item == 'T' or item == 'C':
                totalpyrimidinecount += d[item]
        d1 = {'purine': totalpurinecount, 'pyrimidine': totalpyrimidinecount}
        st.write("Purine and pyrimidine count:", d1)

        # Percentages
        def percentage(dna_seq):
            d2 = {}
            totalcodons = len(dna_seq)
            if totalcodons == 0:
                d2['purinepercentage'] = 0
                d2['pyrimidinepercentage'] = 0
                return d2
            totalpurinecount = dna_seq.count('A') + dna_seq.count('G')
            totalpyrimidinecount = dna_seq.count('T') + dna_seq.count('C')
            d2['purinepercentage'] = (totalpurinecount / totalcodons) * 100
            d2['pyrimidinepercentage'] = (totalpyrimidinecount /
                                          totalcodons) * 100
            return d2

        percentages = percentage(dna)
        st.write("Percentages:", percentages)

        # Reverse Complement
        st.subheader("Reverse Complement")
        dna_list = list(dna)
        new_list = []
        for base in dna_list:
            if base == 'A':
                new_list.append('T')
            elif base == 'T':
                new_list.append('A')
            elif base == 'C':
                new_list.append('G')
            elif base == 'G':
                new_list.append('C')
        new_list.reverse()
        reverse_complement = "".join(new_list)
        st.write("Reverse complement:", reverse_complement)

        # Palindrome Check
        def palindrome_check(sequence):
            return sequence == sequence[::-1]

        st.subheader("Restriction Site Analysis")

        def find_restriction_sites(dna, min_length=4, max_length=8):
            sites = []
            for i in range(min_length, max_length + 1):
                for j in range(len(dna) - i + 1):
                    sub_seq = dna[j:j + i]
                    if palindrome_check(sub_seq):
                        sites.append((j, sub_seq))
            return sites

        restriction_sites = find_restriction_sites(dna)
        if restriction_sites:
            st.write("Found restriction sites:")
            for pos, site in restriction_sites:
                st.write(f"Position {pos}: {site}")
        else:
            st.write("No restriction sites found")

        # GC Content Analysis
        st.subheader("GC Content Analysis")

        def gc_content(dna, window_size):
            dna = dna.upper()
            x, y = [], []
            for i in range(len(dna) - window_size + 1):
                sub_seq = dna[i:i + window_size]
                gc_count = sub_seq.count('G') + sub_seq.count('C')
                gc_percentage = (gc_count / window_size) * 100
                x.append(i)
                y.append(gc_percentage)
            return x, y

        window_size = st.slider("Select window size for GC content analysis:",
                                3, 10, 5)
        x, y = gc_content(dna, window_size)

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(x, y, marker='o')
        ax.set_title("GC content across the DNA sequence")
        ax.set_xlabel("Window position")
        ax.set_ylabel("GC percentage")
        st.pyplot(fig)

        # SNP Detection Section
        st.subheader("SNP Detection")
        reference_file = st.file_uploader("Upload Reference Sequence",
                                          type=["txt", "fasta"],
                                          key="ref")
        query_file = st.file_uploader("Upload Query Sequence",
                                      type=["txt", "fasta"],
                                      key="query")

        def snp_detector(referencefile, queryfile):
            try:
                reference = read_fatsa(referencefile)
                query = read_fatsa(queryfile)

                if len(reference) != len(query):
                    st.warning("Warning: The sequences are not of same length")
                elif len(reference) == 0 or len(query) == 0:
                    st.error("Both sequences are empty")
                    return []

                snp_positions = []
                min_length = min(len(reference), len(query))
                for i in range(min_length):
                    if reference[i] != query[i]:
                        snp_positions.append(i)
                return snp_positions

            except Exception as e:
                st.error(f"Error reading files: {e}")
                return []

        if reference_file and query_file:
            snps = snp_detector(reference_file, query_file)
            if snps:
                st.write(f"SNPs found at positions: {snps}")
            else:
                st.write("No SNPs detected.")

        # Transcription
        st.subheader("DNA to mRNA Transcription")

        def transcribe(dna):
            return dna.replace('T', 'U')

        mrna = transcribe(dna)
        st.write("mRNA sequence:", mrna)

    except FileNotFoundError:
        st.error("File not found. Please check and try again.")
    except Exception as e:
        st.error(f"Error reading file: {e}")

else:
    st.info("Please upload a FASTA file to begin DNA analysis.")
