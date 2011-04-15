
task :env do
  require 'bio'
  TMP = "tmp"
end

task :tmp_dir => :env do
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
