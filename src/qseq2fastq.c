/*
 * qseq2fastq converts Illumina qseq files to fastq format.
 * Written by Andreas Hauser <hauser@genzentrum.lmu.de>, <andy@splashground.de>.
 * Version 1.0 was based on qseq2fastq.pl from http://qseq2fastq.sourceforge.net/.
 * License: GPL Version 3, http://www.gnu.org/licenses/gpl-3.0-standalone.html
 *
 */

#define _GNU_SOURCE 1

#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#define MAX_LINE_LENGTH 4096
#define MAX_BARCODES 1024

int qseq2fastq_demultiplex(FILE *input_qseq_file, FILE *demultiplex_reads_file, char **barcodes, int barcodes_num, FILE** barcode_output_files, int all)
{
  if(input_qseq_file && demultiplex_reads_file)
  {
    char line[MAX_LINE_LENGTH], run[MAX_LINE_LENGTH], read[MAX_LINE_LENGTH], quality_string[MAX_LINE_LENGTH], qfilterc;
    int n, e[12], qfilter;
    int barcode_len[barcodes_num];
    for(int i = 0; i < barcodes_num; i++)
      barcode_len[i] = strlen(barcodes[i]);
    while((n = fscanf(demultiplex_reads_file, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d",
                        run, &e[1], &e[2], &e[3], &e[4], &e[5], &e[6], &e[7], read, quality_string, &qfilter)) > 0)
    {
      if(n != 11)
      {
        fprintf(stderr, "wrong numbers of elements in barcode file line\n");
        return 1;
      }

      if(all || qfilter == 1)
      {
        /* Convert quality filter from 0/1  to N/Y */
        qfilterc = 'N';
        if(qfilter == 1)
          qfilterc = 'Y';
      }
      else
      {
        n = fscanf(input_qseq_file, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d", run, &e[1], &e[2], &e[3], &e[4], &e[5], &e[6], &e[7], read, quality_string, &qfilter);
        continue;
      }

      FILE *output_file = stdout;
      for(int i = 0; i < barcodes_num; i++)
        if(strncmp(barcodes[i], read, barcode_len[i]) == 0)
          output_file = barcode_output_files[i];

      n = fscanf(input_qseq_file, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d", run, &e[1], &e[2], &e[3], &e[4], &e[5], &e[6], &e[7], read, quality_string, &qfilter);
      if(n != 11)
      {
        fprintf(stderr, "wrong numbers of elements in qseq file line\n");
        return 1;
      }

      /* Convert from Illumina encoding of Phred to Sanger encoding */
      for(int c = 0; quality_string[c] != '\0'; c++)
        quality_string[c] = quality_string[c] - 64 + 33;

      /* Convert unknown base from . to N */
      for(int c = 0; read[c] != '\0'; c++)
        if(read[c] == '.')
          read[c] = 'N';

      fprintf(output_file, "@%d:%d:%d:%d:%c\n%s\n+%d:%d:%d:%d:%c\n%s\n", e[2], e[3], e[4], e[5], qfilterc, read, e[2], e[3], e[4], e[5], qfilterc, quality_string);
    }
    return 0;
  }
  return 1;
}

int qseq2fastq(FILE *fp, int all)
{
  if(fp)
  {
    char line[MAX_LINE_LENGTH], run[MAX_LINE_LENGTH], read[MAX_LINE_LENGTH], quality_string[MAX_LINE_LENGTH], qfilterc;
    int n, e[12], qfilter;
    while((n = fscanf(fp, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d",
                        run, &e[1], &e[2], &e[3], &e[4], &e[5], &e[6], &e[7], read, quality_string, &qfilter)) > 0)
    {
      if(n != 11)
      {
        fprintf(stderr, "wrong numbers of elements in line");
        return 1;
      }

      if(all || qfilter == 1)
      {
        /* Convert quality filter from 0/1  to N/Y */
        qfilterc = 'N';
        if(qfilter == 1)
          qfilterc = 'Y';
      }
      else
        continue;

      /* Convert from Illumina encoding of Phred to Sanger encoding */
      for(int c = 0; quality_string[c] != '\0'; c++)
        quality_string[c] = quality_string[c] - 64 + 33;

      /* Convert unknown base from . to N */
      for(int c = 0; read[c] != '\0'; c++)
        if(read[c] == '.')
          read[c] = 'N';

      printf("@%d:%d:%d:%d:%c\n%s\n+%d:%d:%d:%d:%c\n%s\n", e[2], e[3], e[4], e[5], qfilterc, read, e[2], e[3], e[4], e[5], qfilterc, quality_string);
    }
    return 0;
  }
  return 1;
}


void usage()
{
  puts("USAGE: qseq2fastq [-a] [-d BARCODE_QSEQ_FILE -b BARCODE[,BARCODE]* -o OUTPUT_BASE] INPUT_QSEQ_FILES ... > OUTPUT_FASTQ_FILE");
  puts("Converts Illumina qseq files to fastq format. Use - for STDIN as input.");
  puts("\t-a\t\t\tconvert all (don't filter reads flagged with 0/N)");
  puts("\t-d BARCODE_QSEQ_FILE\tdemultiplex according to the qseq file containing the barcode reads");
  puts("\t-b BARCODE[,BARCODE]*\twriting to one output file per barcode and STDOUT for the rest.\n"
                        "\t\t\t\tNo spaces! Barcode matches when it is a valid prefix of the barcode read.");
  puts("\t-o OUTPUT_BASE\t\tprefix used for constructing output filenames");
  puts("\t-h\t\t\tshows this message");
  puts("\nReport bugs to Andreas hauser <andy@splashground.de>.");
}


int main(int argn, char **argv)
{
  char *filename;
  FILE *fp;
  int err, c, all_flag = 0;
  char* barcodes[MAX_BARCODES];
  int barcodes_num = 0;
  char* demultiplex_filename;
  char* output_base = "";
  char* barcode_string;

  while ((c = getopt (argn, argv, "ab:d:ho:")) != -1)
    switch (c)
    {
      case 'a':
        all_flag = 1;
        break;
      case 'd':
        demultiplex_filename = optarg;
        break;
      case 'b':
        barcode_string = strdup(optarg);
        for (int i = 0;  ; i++, barcode_string = NULL)
        {
          char *saveptr1;
          char* token = strtok_r(barcode_string, ",", &saveptr1);
          if(token != NULL)
            barcodes[barcodes_num++] = token;
          else
            break;
        }
        break;
      case 'o':
        output_base = optarg;
        break;
      case 'h':
        usage();
        break;
      case '?':
        if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        return 1;
    }

  if(optind == argn)
  {
    puts("ERROR: Need input file (or \"-\" for STDIN)!");
    usage();
    return 1;
  }

  if(barcodes_num > 0 && demultiplex_filename == NULL)
  {
    puts("ERROR: Barcodes (-b) need demultiplex file (-d)!");
    usage();
    return 1;
  }

  for(int i = optind; i < argn; i++)
  {
    filename = argv[i];
    if(filename[0] == '-' && filename[1] == '\0')
      fp = stdin;
    else
      fp = fopen(filename, "r");

    if(barcodes_num > 0)
    {
      FILE* barcode_output_files[barcodes_num];
      FILE* demultiplex_reads_file = fopen(demultiplex_filename, "r");
      if(demultiplex_reads_file == NULL) { perror("Couldn't open demultiplex file."); return 1; }
      char filename[4096];
      for(int i = 0; i < barcodes_num; i++)
      {
        filename[0] = '\0';
        strcat(filename, output_base);
        strcat(filename, barcodes[i]);
        strcat(filename, ".fastq");
        barcode_output_files[i] = fopen(filename, "a");
        if(barcode_output_files[i] == NULL) { puts(filename); perror("Couldn't open barcode output file."); return 1; }
      }
      qseq2fastq_demultiplex(fp, demultiplex_reads_file, barcodes, barcodes_num, barcode_output_files, all_flag);
    }
    else
      err = qseq2fastq(fp, all_flag);

    if(err > 0)
      fprintf(stderr, " error in filename: %s\n", filename);
  }
  return err;
}


