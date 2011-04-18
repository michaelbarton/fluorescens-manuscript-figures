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

desc "Use phmmer to search all genes against each other"
task :hmmer => [:env,:gene_database] do
  Dir.chdir(TMP) do
    `phmmer --noali --cpu 4 --tblout hits.tab genes.faa genes.faa > phmmer.txt`
  end
end

# "Private" tasks

task :gene_database => :env do
  database = Array.new

  File.open(File.join(TMP,"genes.faa"),"w") do |out|
    Dir['data/reference/gene/**/*.fna'].each do |file|
    source = file.split('/')[-3,3].join('_').gsub(".fna","")
      Bio::FlatFile.auto(file).each do |gene|
        protein = Bio::Sequence::NA.new(gene.seq).translate
        out.puts protein.to_fasta(database.count)

        database << source
      end
    end
  end
  File.open(File.join(TMP,"source.yml"),"w") do |out|
    out.print YAML.dump(database)
  end
end

task :env do
  require 'bio'
  require 'yaml'
  TMP = "tmp"
  FileUtils.rm_rf TMP if File.exists? TMP
  FileUtils.mkdir TMP
end

task :scaffold do
  Dir.chdir('data/genome/assembly') do
    `scaffolder sequence genome.scaffold.yml draft.fna > assembly.fna`
  end
end

task :tmp_dir => [:env,:scaffold] do
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
