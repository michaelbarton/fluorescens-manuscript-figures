task :env do
  require 'pp'

  @tmp    = "tmp"
  @genome =  "genome.fna"
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
  #task :size => [:env,:tmp,'fasta:all'] do
  task :size => [:env] do
    require 'bio'

    File.open('data/genome_size.csv','w') do |out|
      out.puts %W|species source size| * ','
      Dir['tmp/*.fna'].each do |file|
        dna = Bio::FlatFile.auto(file).first.to_biosequence
        source = dna.definition =~ /genome/ ? 'genome' : 'plasmid'
        out.puts([dna.definition.split('_').first,source,dna.seq.length] * ',')
      end
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
