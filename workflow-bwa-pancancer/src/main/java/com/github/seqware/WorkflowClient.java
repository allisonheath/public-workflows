package com.github.seqware;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;

import java.util.ArrayList;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class WorkflowClient extends OicrWorkflow {

    // GENERAL
    // comma-seperated for multiple bam inputs
    String inputBamPaths = null;
    ArrayList<String> bamPaths = new ArrayList<String>();
    // used to download with gtdownload
    String gnosInputFileURL = null;
    String gnosUploadFileURL = null;
    String gnosKey = null;
    // number of splits for bam files, default 1=no split
    int bamSplits = 1;
    String reference_path = null;
    String outputPrefix = null;
    String outputDir = null;
    String dataDir = "data/";
    String outputFileName = "merged_output.bam";
    //BWA
    String bwaAlignMemG = "8";
    String bwaSampeMemG = "8";
    String RGID;
    String RGLB;
    String RGPL;
    String RGPU;
    String RGSM;
    String additionalPicardParams;
    String picardReadGrpMem = "8";
    int readTrimming; //aln
    int numOfThreads; //aln 
    int pairingAccuracy; //aln
    int maxInsertSize; //sampe
    String readGroup;//sampe
    String bwa_aln_params;
    String bwa_sampe_params;
    
    private String picardMaxHeap;
    private String maxMemory;
    private String samtoolsVersion;
    private String picardToolsVersion;
    private String bwaNumberOfThreads;
    private String referencePath;
    
    private static final String LOGS_DIRECTORY = "logs/";
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        /*
        This workflow isn't using file provisioning since 
        we're using GeneTorrent. So this method is just being
        used to setup various variables.
        */
        
        try {
            
            inputBamPaths = getProperty("input_bam_paths");
            for(String path : inputBamPaths.split(",")) {
                bamPaths.add(path);
            }
            gnosInputFileURL = getProperty("gnos_input_file_url");
            gnosUploadFileURL = getProperty("gnos_output_file_url");
            gnosKey = getProperty("gnos_key");
            reference_path = getProperty("input_reference");
            outputDir = this.getMetadata_output_dir();
            outputPrefix = this.getMetadata_output_file_prefix();
            RGID = getProperty("RGID");
            RGLB = getProperty("RGLB");
            RGPL = getProperty("RGPL");
            RGPU = getProperty("RGPU");
            RGSM = getProperty("RGSM");
            bwaAlignMemG = getProperty("bwaAlignMemG") == null ? "8" : getProperty("bwaAlignMemG");
            bwaSampeMemG = getProperty("bwaSampeMemG") == null ? "8" : getProperty("bwaSampeMemG");
            picardReadGrpMem = getProperty("picardReadGrpMem") == null ? "8" : getProperty("picardReadGrpMem");
            additionalPicardParams = getProperty("additionalPicardParams");
            
            setPicardMaxHeap(getProperty("picard_max_heap"));
            setMaxMemory(getProperty("max_memory"));
            setSamtoolsVersion(getProperty("samtools_version"));
            setPicardToolsVersion(getProperty("picard_tools_version"));
            setBwaNumberOfThreads(getProperty("bwa_num_threads"));
            setReferencePath(getProperty("input_reference"));


        } catch (Exception e) {
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.SEVERE, null, e);
            System.exit(1);
        }

        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        // creates the final output 
        this.addDirectory(dataDir);
        
        this.addDirectory(LOGS_DIRECTORY);
    }

    @Override
    public void buildWorkflow() {
        
        ArrayList<Job> bamJobs = new ArrayList<Job>();
        
        // DOWNLOAD DATA
        // let's start by downloading the input BAMs
        Job gtDownloadJob = this.getWorkflow().createBashJob("gtdownload");
        gtDownloadJob.getCommand().addArgument("gtdownload")
                .addArgument("-c "+gnosKey)
                .addArgument("-v -d")
                .addArgument(gnosInputFileURL);
        
        // TODO for loop here
        int numBamFiles = bamPaths.size();
        for (int i=0; i<numBamFiles; i++) {
            
            String file = bamPaths.get(i);
            
            /* samtools view -u -h -f 1 ~/data/SRR062634.bam 2> /tmp/samtools_log.err | 
             * java -Xmx2G -jar SamToFastq.jar INPUT=/dev/stdin INTERLEAVE=true TMP_DIR=/tmp/ FASTQ=/dev/stdout QUIET=true 2> /tmp/samtofastq_log.err | 
             * /usr/bin/bwa mem -p -M -T 0 -R '@RG\tID:SRR062634\tLB:2845856850\tSM:HG00096\tPI:206\tCN:WUGSC\tPL:ILLUMINA\tDS:SRP001294' -t 4 /Users/siakhnin/data/reference/genome.fa.gz 2> /tmp/bwa_log.err - | 
             * java -Xmx2G -jar FixMateInformation.jar INPUT=/dev/stdin OUTPUT=sorted_test.bam SORT_ORDER=coordinate QUIET=true 2> /tmp/fixmate_log.err */
        
            String baseDir = this.getWorkflowBaseDir();
            String samPath = baseDir + "/bin/samtools-" + getProperty("samtools-version") + "/";
            String picardPath = baseDir + "/bin/picard-tools-" + getProperty("picard-tools-version") + "/";
            
            
            Job myValidationJob = this.getWorkflow().createBashJob("bam_validation_" + i);
            myValidationJob.addParent(gtDownloadJob);
            myValidationJob.setMaxMemory(getMaxMemory());
            
            Command validationJobCommand = myValidationJob.getCommand();
            validationJobCommand.addArgument("java -Xmx" + getPicardMaxHeap() + " -jar " + picardPath + "ValidateSamFile.jar INPUT=" + file + " OUTPUT=" + LOGS_DIRECTORY + file + ".validation.log");
            
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.INFO, null, validationJobCommand.toString());
            
            
            Job myBwaJob = this.getWorkflow().createBashJob("bwa_mem_" + i);
            myBwaJob.addParent(myValidationJob);
            myBwaJob.setMaxMemory(getMaxMemory());
            
            Command bwaJobCommand = myBwaJob.getCommand();
            
            
            //Samtools filters out unpaired reads -u = uncompressed, -h = with header, -f 1 = flag for paired reads
            bwaJobCommand.addArgument(samPath + "samtools view ");
            bwaJobCommand.addArgument(" -u -h -f 1 " + file + " 2> ./" + file + ".samtools.err ");
            
            //Picard SamToFastq
            bwaJobCommand.addArgument("| java -Xmx" + getPicardMaxHeap() + " -jar " + picardPath + "SamToFastq.jar");
            bwaJobCommand.addArgument("INPUT=/dev/stdin INTERLEAVE=true TMP_DIR=./ FASTQ=/dev/stdout QUIET=true 2> " + LOGS_DIRECTORY + file + ".samtofastq.err ");
            
            //BWA mem
            bwaJobCommand.addArgument("| " + baseDir + "/bin/bwa mem ");
            bwaJobCommand.addArgument("-p -M -T 0 -t " + getProperty("bwa_num_threads") + " " + getProperty("input_reference") + " 2> ." + LOGS_DIRECTORY + file + ".bwa.err - ");
            
            //Fix Mate Info and sort by coordinate
            bwaJobCommand.addArgument("java -Xmx" + getPicardMaxHeap() + " -jar " + picardPath + "FixMateInformation.jar ");
            bwaJobCommand.addArgument("INPUT=/dev/stdin OUTPUT=" + dataDir + file + ".pcap.bam SORT_ORDER=coordinate QUIET=true 2> " + LOGS_DIRECTORY + file + ".fixmate.err ");
            
            
            Logger.getLogger(WorkflowClient.class.getName()).log(Level.INFO, null, bwaJobCommand.toString());
            
            Job myUploadJob = this.getWorkflow().createBashJob("upload");
            Command myUploadCommand = myUploadJob.getCommand();
            myUploadCommand.addArgument("perl " + this.getWorkflowBaseDir() + "/scripts/upload_data.pl");
            myUploadCommand.addArgument("--bam " + dataDir + file + ".pcap.bam");
            myUploadCommand.addArgument("--key " + gnosKey);
            myUploadCommand.addArgument("--url " + gnosUploadFileURL);
            myUploadJob.addParent(myBwaJob);
                      
        }
    }

    public String parameters(final String setup) {

        String paramCommand = null;
        StringBuilder a = new StringBuilder();

        try {
            if (setup.equals("aln")) {

                if (!getProperty("numOfThreads").isEmpty()) {
                    numOfThreads = Integer.parseInt(getProperty("numOfThreads"));
                    a.append(" -t ");
                    a.append(numOfThreads);
                    a.append(" ");
                }

                if (!getProperty("bwa_aln_params").isEmpty()) {
                    bwa_aln_params = getProperty("bwa_aln_params");
                    a.append(" ");
                    a.append(bwa_aln_params);
                    a.append(" ");
                }
                paramCommand = a.toString();
                return paramCommand;
            }

            if (setup.equals("sampe")) {

                if (!getProperty("maxInsertSize").isEmpty()) {
                    maxInsertSize = Integer.parseInt(getProperty("maxInsertSize"));
                    a.append(" -a ");
                    a.append(maxInsertSize);
                    a.append(" ");
                }

                if (!getProperty("readGroup").isEmpty()) {
                    a.append(" -r ");
                    a.append(readGroup);
                    a.append(" ");
                }

                if (!getProperty("bwa_sampe_params").isEmpty()) {
                    bwa_sampe_params = getProperty("bwa_sampe_params");
                    a.append(" ");
                    a.append(bwa_sampe_params);
                    a.append(" ");
                }
                paramCommand = a.toString();
                return paramCommand;
            }

        } catch (Exception e) {
        }
        return paramCommand;
    }

	public String getPicardMaxHeap() {
		return picardMaxHeap;
	}

	public void setPicardMaxHeap(String picardMaxHeap) {
		this.picardMaxHeap = picardMaxHeap;
	}

	public String getMaxMemory() {
		return maxMemory;
	}

	public void setMaxMemory(String maxMemory) {
		this.maxMemory = maxMemory;
	}

	public String getSamtoolsVersion() {
		return samtoolsVersion;
	}

	public void setSamtoolsVersion(String samtoolsVersion) {
		this.samtoolsVersion = samtoolsVersion;
	}

	public String getPicardToolsVersion() {
		return picardToolsVersion;
	}

	public void setPicardToolsVersion(String picardToolsVersion) {
		this.picardToolsVersion = picardToolsVersion;
	}

	public String getBwaNumberOfThreads() {
		return bwaNumberOfThreads;
	}

	public void setBwaNumberOfThreads(String bwaNumberOfThreads) {
		this.bwaNumberOfThreads = bwaNumberOfThreads;
	}

	public String getReferencePath() {
		return referencePath;
	}

	public void setReferencePath(String referencePath) {
		this.referencePath = referencePath;
	}
}
