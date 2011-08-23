task :env do
  require 'fastercsv'
  require 'pp'
  require 'set'
  require 'bio'

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

  desc "Perform all-to-all hmmer gene search"
  task :search => [:env,'fasta:database'] do
    Dir.chdir(@tmp) do
      `phmmer --noali --cpu 7 -E #{0.001} --tblout hits.tab genes.faa genes.faa > phmmer.txt`
    end

    hits = Hash.new{|h,k| h[k] = []}
    File.open('tmp/hits.tab').reject{|i| i[0,1] == '#'}.map do |i|
      entries = i.split
      match, query, score = entries[0], entries[2], entries[4].to_f

      hits[query] << {:match => match, :score => score}
    end

    delim = "\t"
    File.open('data/search.tab','w') do |out|
      out.puts %w|query match score| * delim
      hits.keys.sort.each do |query|
        hits[query].sort_by{|i| i[:score]}.each do |hit|
          out.puts [query, hit[:match], hit[:score]] * delim
        end
      end
    end
  end

end

namespace :orthologs do

  desc "Calculate ortholog clusters"
  task :all => [:tmp,'data:search',:cluster,:parse]

  task :cluster => [:env] do
    Dir.chdir(@tmp) do
      `mcxload -abc ../data/search.tab -o graph.txt -write-tab labels.txt`
      `mcl graph.txt -o clusters.txt -use-tab labels.txt --force-connected=y`
      `clmformat -icl clusters.txt -imx graph.txt -dir . -dump cluster-scores.txt --dump-measures`
      FileUtils.cp "clusters.txt","../data/orthologs/"
    end
  end

  task :parse => [:env] do
    clusters = File.open(@tmp + '/clusters.txt').map{ |line| line.split }
    species = YAML.load(File.read('data/species.yml'))

    (species.values.flatten - clusters.flatten).each do |id|
      clusters << [id]
    end

    File.open('data/orthologs/clusters.yml','w') do |out|
      out.puts YAML.dump(clusters)
    end

    invert = Hash.new
    species.each do |specie,genes|
      genes.each do |gene|
        invert[gene] = specie
      end
    end

    species_clusters = clusters.map do |cluster|
      cluster.map{|gene| invert[gene] }.uniq
    end

    File.open('data/orthologs/species_clusters.yml','w') do |out|
      out.puts YAML.dump(species_clusters)
    end

  end

end

namespace :pili do

  task :all => [:tmp,'fasta:genes']

  task :search => :env do
    search  = Set.new # List of gene IDs to search
    expired = Set.new # Already searched genes

    results = Array.new

    # Seed search set with initial data
    search += YAML.load(File.read('data/pili/plasmid_pili_gene_ids.yml')).keys

    until(search.empty?) do
      puts "Searching #{format("%3d",search.size)} genes"
      hits = phmmer(search)

      results += hits

      expired |= search
      search.clear
      search  |= (hits.map{|hit| hit[:match]} - expired.to_a)
    end
    File.open('data/pili/adjacency.csv','w') do |out|
      out.puts %w|query match score| * ','
      results.sort_by{|row| [row[:query],row[:match]]}.each do |row|
        out.puts [row[:query],row[:match],row[:score]] * ','
      end
    end
  end

  task :create_db => :env do
    database do |db|
      File.open('data/search.tab').each_with_index do |line,i|
        query, match, score = line.split
        db[i] = {:query => query, :match => match, :score => score.to_f}
      end
    end
  end

  task :plot_data => :env do
    names = YAML.load(File.read("data/pili/plasmid_pili_gene_ids.yml"))
    File.open('data/pili/plot_data.csv','w') do |out|
      out.puts %w|gene,x,y,name,source| * ','
      FCSV.foreach("data/pili/scaling.csv",:headers => true) do |row|

        # Label with pili gene name
        row << (names['plasmid'][row[0]] || names['genome'][row[0]] || "")

        # Label with source
        if    names['plasmid'][row[0]]
          row << 'plasmid'
        elsif names['genome'][row[0]]
          row << 'genome'
        else
          row << ''
        end

        out.puts row
      end
    end
  end

  def sequences(ids)
    Bio::FlatFile.auto("tmp/genes.faa").select do |entry|
      ids.include? entry.definition.split.first
    end
  end

  def database(&block)
    require 'rufus/tokyo'
    file = 'tmp/search.tct'
    db = Rufus::Tokyo::Table.new(file)
    yield db
    db.close
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

  task :database => [:env] do
    require 'bio'
    database = Hash.new{|h,k| h[k] = [] }

    File.open(File.join(@tmp,"genes.faa"),"w") do |out|

      Dir['data/reference/gene/**/*.fna'].each do |file|
      source = file.split('/')[-3,3].join('_').gsub(".fna","")
        Bio::FlatFile.auto(file).each do |gene|
          name = genebank_id(gene.definition.split.first)
          protein = Bio::Sequence::NA.new(gene.seq).translate
          out.puts protein.to_fasta(name)
          database[source] << name
        end
      end

      source = "R124"
      Bio::FlatFile.auto("data/genome/annotation/genes.fna").each do |gene|
        name = gene.definition.split.first
        protein = Bio::Sequence::NA.new(gene.seq).translate
        out.puts protein.to_fasta(name)
        database[source] << name
      end

    end
    File.open(File.join("data/species.yml"),"w") do |out|
      out.print YAML.dump(database)
    end
  end

  def genebank_id(string)
    string =~ /([^_]+)\.\d$/
    Regexp.last_match(1)
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
