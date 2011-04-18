desc "Calculates ortholog clusters"
task :orthologs => [:database,:hmmer,:parse].map{|i| "orthologs:{i}"}
namespace :orthologs do

  task :cluster => :env do
    Dir.chdir(TMP) do
      `mcxload -abc network.csv -o graph.txt -write-tab labels.txt`
      `mcl graph.txt -o clusters.txt -use-tab labels.txt --force-connected=y`
      `clmformat -icl clusters.txt -imx graph.txt -dir . -dump cluster-scores.txt --dump-measures`
    end
  end

  task :parse => :env do
    Dir.chdir(TMP) do
      File.open('network.csv','w') do |out|
        File.open('short.tab','r').each do |line|
          next if line[0,1] == '#'
          tokens = line.split

          a,b = tokens[0],tokens[2]
          next if a == b
          next if tokens[4].to_f > 1e-3 # Test for homology using e-value

          out.puts [a,b].join(' ')
        end
      end
    end
  end

  task :hmmer => :env do
    Dir.chdir(TMP) do
      `phmmer --noali --cpu 4 --tblout hits.tab genes.faa genes.faa > phmmer.txt`
    end
  end

  task :database => :env do
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
end

desc "Run genome alignments"
task :alignments => [:cp_genomes,:nucmer].map{|i| "alignments:{i}"}
namespace :alignments do

  task :nucmer => [:env,:scaffold] do
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

  task :cp_genomes => :tmp_dir do
    Dir['data/reference/genomes/**/*.gb'].each do |f|
      dna = Bio::FlatFile.auto(f).first.to_biosequence

      target = f.split('/')[-3,3] * '_'
      target.gsub!(".gb",".fna")
      File.open(File.join(TMP,target),"w") do |out|
        out.print dna.output(:fasta,:header => target.gsub(".fna",""))
      end
    end
  end

end

# "Private" tasks

task :env do
  require 'bio'
  require 'yaml'
  TMP = "tmp"
end

task :scaffold do
  Dir.chdir('data/genome/assembly') do
    `scaffolder sequence genome.scaffold.yml draft.fna > assembly.fna`
  end
end

task :tmp_dir => :env do
  FileUtils.rm_rf TMP if File.exists? TMP
  FileUtils.mkdir TMP
end
