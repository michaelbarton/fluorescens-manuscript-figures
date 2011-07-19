task :env do
  require 'pp'

  @tmp    = "tmp"
  @genome =  "fluorescens_r124_genome.fna"
end

task :tmp  => :env do
  @tmp = File.expand_path File.join(File.dirname(__FILE__),'tmp')
  FileUtils.rm_rf @tmp if File.exists? @tmp
  FileUtils.mkdir @tmp
end

namespace :data do

  desc "Align genomes to the scaffold"
  task :synteny => [:env,:tmp,'fasta:all'] do
    Dir.chdir(@tmp) do
      pids = Dir['*_*.fna'].map do |r|
        ref = r.gsub('.fna','')
        pid = fork{
        `nucmer --prefix=#{ref} --maxmatch #{ref}.fna #{@genome} 2> /dev/null`
        `show-coords -THr #{ref}.delta > #{ref}.coords`
        }
        pid
      end
    # Wait for all jobs to finish
    pids.each {|p| Process.waitpid p}
    end
    `sort #{@tmp}/*.coords -o data/alignment/nucmer.coords`
  end

  desc "Calculate genome sizes"
  task :size => [:env,:tmp,'fasta:scaffold'] do
    require 'bio'
    require 'fastercsv'

    files = Dir['data/reference/genomes/**/*.gb']
    r124_genes = "data/genome/annotation/gene_list.csv"

    data = files.map do |file|
      datum = file.gsub('.gb','').split('/')[-3,3]
      entry = Bio::FlatFile.auto(file).first
      datum << entry.seq.length

      count = 0
      entry.each_gene{ count += 1}
      datum << count
    end

    r124 = %w|fluorescens R124 genome|
    r124 << Bio::FlatFile.auto('tmp/fluorescens_r124_genome.fna').first.seq.length
    r124 << FCSV.open(r124_genes,:headers => true).inject(0) do |count,row|
      count += 1 if row['Scaffold Name'] =~ /scaffold|contig00001/
      count
    end
    data << r124

    FasterCSV.open('data/genome_size.csv','w') do |out|
      out << %w|Pseudomonas strain source genome_size gene_count|
      data.each{|row| out << row}
    end
  end

end

namespace :fasta do

  task :all => [:reference,:scaffold]

  task :reference => [:env,:tmp] do
    require 'bio'

    files = Dir['data/reference/genomes/**/*.gb']
    files.each do |f|
      dna = Bio::FlatFile.auto(f).first.to_biosequence

      target = f.split('/')[-3,3] * '_'
      target.gsub!(".gb",".fna")
      File.open(File.join(@tmp,target),"w") do |out|
        out.print dna.output(:fasta,:header => target.gsub(".fna",""))
      end
    end
  end

  task :scaffold => [:env,:tmp] do
    Dir.chdir('data/genome/assembly') do
      `scaffolder sequence genome.scaffold.yml draft.fna --definition="fluorescens_r124_genome" > #{@tmp}/#{@genome}`
    end
  end

end

namespace :plot do

  task :all => [:genome_size,:genome_dot_plot]

  task :genome_size => 'data:size' do
    `Rscript plot/genome_size_vs_annotations.R > /dev/null`
  end

  task :genome_dot_plot => 'data:synteny' do
    `Rscript plot/genome_dot_plot.R > /dev/null`
  end

end
