#! /usr/bin/env nextflow

// Copyright (C) 2010 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.mem           = 8
params.cpu           = 4
params.output_folder = "bam2peaks"
params.input_file    = null
params.ref           = "hg38"
params.gencode       = null
params.ignoreDuplicates = null
params.extsize = 150
params.broad = null
params.mode = "atac"
params.saturation = null

params.help = null

log.info "" 
log.info "--------------------------------------------------------"
log.info "  mutect-nf 2.2.0: Mutect pipeline for somatic variant calling with Nextflow DSL2"
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""


if (params.help) {
log.info '-------------------------------------------------------------'
    log.info ' USAGE  '
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'nextflow run mutect.nf --tumor_bam_folder tumor_BAM/ --normal_bam_folder normal_BAM/ --ref ref.fasta [OPTIONS] '
    log.info ''
    log.info 'Mandatory arguments:'
    log.info ''
    log.info '    --input_file         FILE                    input tabulation-separated values file with columns sample (sample name),'
    log.info '                                                 tag (short name for figures), bam (bam file path) and group (group)'
    log.info ''
    log.info '    --ref                STRING                  Reference fasta file. hg19, hg38 or mm10'
    log.info ''
    log.info 'Optional arguments:'
    log.info ''
    log.info '    --mode               STRING                  There is two mode : atac or chip (default = atac), chip require input sample(s)'    
    log.info '    --cpu                INTEGER                 Number of cpu used (default: 16).'
    log.info '    --mem                INTEGER                 Java memory passed to mutect in GB (default: 16).'
    log.info '    --extsize            INTEGER                 MACS extsize : extendsize of peaks to to fix-sized fragments (default: 150).'    
    log.info '    --output_folder      FOLDER                  Output folder (default: bam2peaks).'
    log.info ''
    log.info 'Flags:'
    log.info ''
    log.info '    --broad                                     Compute broadpeaks instead of narrowpeaks'
    log.info '    --ignoreDuplicates                          Ignore duplicates reads'
    log.info '    --saturation                                Run saturation process'
    log.info ''
    exit 0
}else{
    /* Software information */
    log.info "mem                    = ${params.mem}"
    log.info "cpu                    = ${params.cpu}"
    log.info "mode                   = ${params.mode}"
    log.info "ignoreDuplicates       = ${params.ignoreDuplicates}"
    log.info "output_folder          = ${params.output_folder}"
    log.info "input_file             = ${params.input_file}"
    log.info "ref                    = ${params.ref}"
    log.info "extsize                = ${params.extsize}"
    log.info "broad                  = ${params.broad}"
    log.info "saturation             = ${params.saturation}"
}


/***************************************************************************************/
/************************ handle global parameters *************************************/
/***************************************************************************************/

//rmd_script = file("${baseDir}/saturation.Rmd")
sampleSheet = file(params.input_file)
//println(sampleSheet)

binSize=10

if(params.ref=="hg38") {
    genomeSize=2913022398
    blackListFile=file("$baseDir/blacklist/hg38-blacklist.v2.bed.gz")
}else if(params.ref=="hg19"){
    genomeSize=2864785220
    blackListFile=file("$baseDir/blacklist/hg19-blacklist.v2.bed.gz")
}else if(params.ref=="mm10"){
    genomeSize=2652783500
    blackListFile=file("$baseDir/blacklist/mm10-blacklist.v2.bed.gz")
}else{
    Error : "param.ref must be one of hg19, hg38, mm10"
}

gencode = file(params.gencode)
broad = params.broad ? "--broad" : ""


/***************************************************************************************/
/************************  Process : bamCoverage  **************************************/
/***************************************************************************************/

process bamCoverage{

    memory params.mem+'GB'
	cpus params.cpu

    publishDir "${params.output_folder}/bw", mode: 'copy', pattern:"*.bw"

    input:
        tuple val(sample), path(bam), path(bai)
        path(blackListFile)

    output:
        path("*.bw")

    shell:
        ignoreDuplicates = params.ignoreDuplicates ? "--ignoreDuplicates" : ""
        """
        bamCoverage -b $bam --effectiveGenomeSize $genomeSize --binSize $binSize -o ${sample}.bw $ignoreDuplicates \
            --normalizeUsing RPGC --blackListFileName $blackListFile -p $task.cpus
        """

    stub:
        ignoreDuplicates = params.ignoreDuplicates ? "--ignoreDuplicates" : ""
        """
        echo bamCoverage -b $bam --effectiveGenomeSize $genomeSize --binSize $binSize -o ${sample}.bw \
            $ignoreDuplicates --normalizeUsing RPGC --blackListFileName $blackListFile -p $task.cpus
        touch ${sample}.bw
        """
}

/***************************************************************************************/
/************************  Process : qc_TSS2kb  ****************************************/
/***************************************************************************************/

process qc_TSS2kb{

    memory 8
	cpus 2

    publishDir "${params.output_folder}/QCs", mode: 'copy'

    input:
        path(bw) // all bw required
        path(gencode)

    output:
        path("*png")

    shell:
        bw_ = bw.join(" ")
        """
        computeMatrix reference-point -S $bw_ -R $gencode -a 2000 -b 2000 -p $task.cpus --outFileName matrix_TSS2kb.gz
        #--samplesLabel WT-S1 WT-S2 Cl11-S3 Cl11-S4 Cl11-S5 BAP1-S15-1 BAP1-S15-2 BAP1-S18-1 BAP1-S18-2
        
        plotHeatmap -m matrix_TSS2kb.gz -o matrix_TSS2kb.png 
        """
    
    stub:
        bw_ = bw.join(" ")
        """
        echo computeMatrix reference-point -S $bw_ -R $gencode -a 2000 -b 2000 -p $task.cpus --outFileName matrix_TSS2kb.gz
        touch matrix_TSS2kb.png
        """
}


/***************************************************************************************/
/************************  Process : fragmentSize  *************************************/
/***************************************************************************************/

process fragmentSize{

    memory 8
	cpus 2

    publishDir "${params.output_folder}/QCs", mode: 'copy'

    input:
        tuple val(tags), path(bams), path(bais) // all element are also tuple containing severals elements

    output:
        path("*png")

    shell:
        bams_=bams.join(" ")
        tags_=tags.join(" ")
        """
        bamPEFragmentSize -hist fragmentSize.png -p $task.cpus -T "Fragment size of PE ATAC-seq data" --maxFragmentLength 1000 -b $bams_ \
        --samplesLabel $tags_
        """

    stub:
        bams_=bams.join(" ")
        tags_=tags.join(" ")
        """
        echo bamPEFragmentSize -hist fragmentSize.png -p $task.cpus -T "Fragment size of PE ATAC-seq data" --maxFragmentLength 1000 \
        -b $bams_ --samplesLabel $tags_
        touch fragmentSize.png
        """
}



/***************************************************************************************/
/************************  Process : multiBamSummary  **********************************/
/***************************************************************************************/

process multiBamSummary{

    memory 8
	cpus 2

    publishDir "${params.output_folder}/QCs", mode: 'copy'

    input:
        tuple val(sample), path(bams), path(bais) // all bams

    output:
        path("*.jpg")

    shell:
        ignoreDuplicates = params.ignoreDuplicates ? "--ignoreDuplicates" : ""
        bams_=bams.join(" ")
        """
        multiBamSummary bins --bamfiles $bams_ $ignoreDuplicates --minMappingQuality 30 -out readCounts.npz
        plotCorrelation --corData readCounts.npz --corMethod spearman --whatToPlot heatmap --plotFile correlationReads.jpg --plotNumbers --removeOutliers
        plotPCA --corData readCounts.npz --plotFile PCAreads.jpg --rowCenter
        """

    stub:
        ignoreDuplicates = params.ignoreDuplicates ? "--ignoreDuplicates" : ""
        bams_=bams.join(" ")
        """
        echo multiBamSummary bins --bamfiles $bams_ $ignoreDuplicates --minMappingQuality 30 -out readCounts.npz
        echo plotCorrelation --corData readCounts.npz --corMethod spearman --whatToPlot heatmap --plotFile correlationReads.jpg --plotNumbers --removeOutliers
        echo plotPCA --corData readCounts.npz --plotFile PCAreads.jpg --rowCenter
        touch correlationReads.jpg PCAreads.jpg
        """
}


/***************************************************************************************/
/************************  Process : fingerPrint  **************************************/
/***************************************************************************************/

process fingerPrint{

    memory params.mem
	cpus params.cpu

    publishDir "${params.output_folder}/QCs", mode: 'copy'

    input:
        tuple val(tags), path(bams), path(bais) // all bams

    output:
        path("*.png")

    shell:
        ignoreDuplicates = params.ignoreDuplicates ? "--ignoreDuplicates" : ""
        bams_=bams.join(" ")
        tags_=tags.join(" ")
        """
        plotFingerprint -b $bams_ --labels $tags --minMappingQuality 30 --skipZeros $ignoreDuplicates -p $task.cpus \
            -T "Fingerprints of different samples" --plotFile fingerprints.png
        """

    stub:
        ignoreDuplicates = params.ignoreDuplicates ? "--ignoreDuplicates" : ""
        bams_=bams.join(" ")
        tags_=tags.join(" ")
        """
        echo plotFingerprint -b $bams_ --labels $tags --minMappingQuality 30 --skipZeros $ignoreDuplicates -p $task.cpus \
            -T "Fingerprints of different samples" --plotFile fingerprints.png
        touch fingerprints.png
        """
}

/***************************************************************************************/
/************************  Process : subSamples  ***************************************/
/***************************************************************************************/

process subSamples{

    memory params.mem
	cpus params.cpu

    publishDir "${params.output_folder}/Counts", mode: 'copy', pattern: "*.count"

    input:
        tuple val(group), val(sample), path(bam), path(bais), path(input_bam), path(input_bai)
        each(percentage)

    output:
        tuple val(group), val(subsample), path("*$subsample*.bam"), path("*$subsample*.bai"), path(input_bam), path(input_bai), emit : bams
        path("*.count"), emit : counts

    shell:
        subsample = sample + "_" + percentage
        """
        if [ ! -e "$input_bam" ]; then touch NO_INPUT1 NO_INPUT2; fi
        samtools view -h -b -s $percentage $bam > ${subsample}.bam
        samtools index ${subsample}.bam
        samtools view -c ${subsample}.bam > ${subsample}.count
        """

    stub:
        subsample = sample + "_" + percentage
        """
        if [ ! -e "$input_bam" ]; then touch NO_INPUT1 NO_INPUT2; fi
        echo samtools view -h -b -s $percentage $bam > ${subsample}.bam
        echo samtools index ${subsample}.bam
        echo samtools view -c ${subsample}.bam > ${subsample}.count
        touch ${subsample}.bam ${subsample}.bai ${subsample}.count
        """
}



/***************************************************************************************/
/************************  Process : macs2  ********************************************/
/***************************************************************************************/

process macs2{

    memory params.mem
	cpus params.cpu

    publishDir "${params.output_folder}/${mode}", mode: 'copy'

    input:
        tuple val(group), val(sample), path(bam), path(bai), path(input_bam), path(input_bai)
        val(mode)

    output:
        tuple val(group), path("*.{narrowPeak,broadPeak}")

    shell:
        keepDuplicates = params.ignoreDuplicates ? "" : "--keep-dup auto"
        inputs_ = (input_bam.name == "NO_INPUT1") ? "" : "-c $input_bam"
        """
        echo $input_bam
        macs2 callpeak -t $bam $inputs_ -n $sample -p 0.01 -g $genomeSize $keepDuplicates $broad --nomodel --extsize $params.extsize
        """

    stub:
        keepDuplicates = params.ignoreDuplicates ? "" : "--keep-dup auto"
        inputs_ = ("$input_bam" == "NO_INPUT1") ? "" : "-c $input_bam"
        """
        echo macs2 callpeak -t $bam $inputs_ -n $sample -p 0.01 -g $genomeSize $keepDuplicates $broad --nomodel --extsize $params.extsize
        touch ${sample}.broadPeak
        """

}


/***************************************************************************************/
/************************  Process : idr  **********************************************/
/***************************************************************************************/


process idr{

    memory params.mem
	cpus params.cpu

    publishDir "${params.output_folder}/Peaks_intersect", mode: 'copy', pattern:"*.bed"

    input:
        tuple val(group), path(peaks)

    output:
        tuple val(group), path("*.bed")

    shell:
        peaks_=peaks.join(" ")
        file_type = peaks_.tokenize('.').last()
        """
        echo $file_type
        idr --samples $peaks_ --input-file-type $file_type -o ${group}.${file_type}
        awk -F "\t" '(\$5 >= 540 )' ${group}.${file_type} | cut -f 1-3 > ${group}.${file_type}_IDR_filter.bed
        """

    stub:
        peaks_=peaks.join(" ")
        """
        touch output_${group}.bed
        """

}

/***************************************************************************************/
/************************  Process : idr  **********************************************/
/***************************************************************************************/

process rmd_report{

	cpus 1

    publishDir "${params.output_folder}", mode: 'copy', pattern:"*html"

    input:
        //path rmd_script
        //path sampleSheet
        path peaks
        path counts

    output:
        path("*.html")

    shell:
        """
        mkdir Saturation_peaks Counts
        mv *_peaks* Saturation_peaks/.
        mv *count Counts/.
        cp ${baseDir}/saturation.Rmd .
        Rscript -e "rmarkdown::render('saturation.Rmd',params=list(sampleSheet = '$sampleSheet'))"
        """

    stub:
        """
        touch saturation.html
        """

}


/***************************************************************************************/
/************************  Workflow : main *********************************************/
/***************************************************************************************/

workflow saturation {

    take: 
        bams4masc2

    main:
        def percentages = Channel.of(0.15,0.3,0.5,0.6,0.7,0.8,0.9,0.95)
        subSamples(bams4masc2,percentages)
        macs2(subSamples.out.bams,"Saturation_peaks")
        foo=macs2.out.map{ it[1] }.toList()
        rmd_report(foo, subSamples.out.counts.collect())

}


workflow {

    /**** INPUTS ****/
    bams = Channel.fromPath(params.input_file).splitCsv(header: true, sep: '\t', strip: true)
        .map{ row ->
        assert(row.sample) : "Error: sample column is missing check your input_file"
        assert(row.bam) : "Error: bam column is missing check your input_file"
        tuple(
            row.group ? row.group : "1", row.input ? row.input : "0", row.sample, row.tag ? row.tag : row.sample,  
            file(row.bam), file("${row.bam}.{bai,crai}")[0]
        )} // group, input, sample, tag, bam, bai


    /**** QCs ****/
    bams4coverage = bams.map{ row-> tuple(row[2], row[4], row[5] ) } // sample, bam, bai
    bams4fragmentSize = bams.map{ row-> tuple(row[3], row[4], row[5] ) }.toList().transpose().toList() // tags, bams, bais

    bamCoverage(bams4coverage,blackListFile)
    
    qc_TSS2kb(bamCoverage.out.collect(), gencode)
    
    fragmentSize(bams4fragmentSize)

    multiBamSummary( bams4coverage.toList().transpose().toList() )

    fingerPrint(bams4fragmentSize)

    /**** MACS2 ****/
    foo = bams.map{ row-> tuple(row[0], row[1], row[2], row[4], row[5] ) } // group, input, sample, bam, bai
        .branch{
            samp: it[1] == "0"
            input: it[1] == "1"
        }

    bams4macs2 = (params.mode=="chip" ? foo.samp.join( foo.samp.combine(foo.input, by:0), remainder: true, by:[0,1,2,3,4]) : foo.samp)
        .map{ row -> tuple( row[0], row[2], row[3], row[4],
                            row[7] ? row[7] : file("NO_INPUT1"), 
                            row[8] ? row[8] : file("NO_INPUT2") 
            )} // group, sample, bam, bai, input_bam, input_bai*/

    macs2(bams4macs2,"Peaks") | groupTuple | idr

    if(params.saturation){
        saturation(bams4macs2)
    }




}


