task :env do
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
      `scaffolder sequence genome.scaffold.yml draft.fna > #{@tmp}/#{@genome}`
    end
  end

end
