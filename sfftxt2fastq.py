#!/usr/bin/env python
# import modules
import argparse, mmap, os, re, sys
from StringIO import StringIO
from Bio import SeqIO

# convert quality scores
def convert_quality_scores(int_scores):
  score_map_dict = {str(i - 33): chr(i) for i in range(33, 73 + 1)}
  str_scores = ""
  for score in int_scores.split():
    str_scores += score_map_dict[score]
  return str_scores

# subroutine to parse sfftxt file
def parse_sfftxt(args):
  # set output filename
  outfile = args.sfftxt.name.rstrip(".sff.txt") + ".fastq"
  open(outfile, "w").close()

  # seek to beginning of file
  mm = mmap.mmap(args.sfftxt.fileno(), 0, prot=mmap.PROT_READ)
  mm.seek(0)

  # iterate over records
  record_offset = None
  num_records = 0
  while not record_offset == -1:
    # find next record
    record_offset = mm.find("\n>") + 1
    if record_offset == -1 or record_offset == 0:
      break
    mm.seek(record_offset)
    lines = []
    # records are 21 lines in length
    for i in range(21):
      lines.append(mm.readline().strip().split("\t"))
      record_dict = {}
      # iterate over lines in record
      for line in lines:
        if len(line) == 1 and line[0] == "":
          pass
        elif len(line) == 1 and not line[0] == "":
          if line[0].startswith(">"):
            record_dict["accession"] = "@" + line[0][1:]
          if re.search(":", line[0]):
            fields = [line.strip() for line in line[0].split(":")]
            record_dict[fields[0]] = fields[1]
        elif len(line) > 2:
          key = line[0][:-1]
          value = " ".join(line[1:])
          record_dict[key] = value
        else:
          record_dict[line[0][:-1]] = line[1]

    # build record
    record = "\n".join([" ".join([record_dict["accession"], record_dict["Run Name"], "length=" + str(len(record_dict["Bases"]))]), record_dict["Bases"].upper(), "+", convert_quality_scores(record_dict["Quality Scores"])])
    record = SeqIO.read(StringIO(record), "fastq")
    num_records += 1 

    # output to file
    SeqIO.write(record, open(outfile, "a"), "fastq")  

  # let user know we're done
  print "Outptted %d sequences to %s." % (num_records, outfile)

# main subroutine
def main():
  parser = argparse.ArgumentParser(prog='sfftxt2fastq', description="Convert from sfftxt format to FASTQ")
  parser.add_argument('sfftxt', type=argparse.FileType('r'), help="input sfftxt file")
  args = parser.parse_args()
  print "Parsing file..."
  parse_sfftxt(args)

# run main if executing from command line
if __name__ == "__main__":
  main()

