desc "Generate genome alignments"
task :alignments => [:env,:tmp_dir] do
  Dir.chdir(TMP) do
    pids = Dir['*.fna'].map do |r|
      ref = r.gsub('.fna','')
      pid = fork{
      `nucmer --prefix=#{ref} --maxmatch #{ref}.fna ../data/genome/assembly/assembly.fna 2> /dev/null && show-coords -THr #{ref}.delta > #{ref}.coords`
      }
      pid
    end
  # Wait for all jobs to finish
  pids.each {|p| Process.waitpid p}
  end
  `cat #{TMP}/*.coords > data/alignment/nucmer.coords`
end

task :env do
  require 'bio'
  TMP = "tmp"
end

task :scaffold do
  Dir.chdir('data/genome/assembly') do
    `scaffolder sequence genome.scaffold.yml draft.fna > assembly.fna`
  end
end

task :tmp_dir => [:env,:scaffold] do
  FileUtils.rm_rf TMP if File.exists? TMP
  FileUtils.mkdir TMP

  files = Dir['data/reference/genomes/**/*.gb']
  files.each do |f|
    dna = Bio::FlatFile.auto(f).first.to_biosequence

    target = f.split('/')[-3,3] * '_'
    target.gsub!(".gb",".fna")
    File.open(File.join(TMP,target),"w") do |out|
      out.print dna.output(:fasta,:header => target.gsub(".fna",""))
    end
  end
end
