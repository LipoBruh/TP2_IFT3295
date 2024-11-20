import * as fs from 'node:fs';
//Saves every score for their corresponding amino acid
//Acts as a dictionnary via the use of the Map class
class AminoScore{
    
    constructor(symbol="*"){
        this.data = this.symbolToAminoAcid(symbol);
        this.map = new Map();
    }

    symbolToAminoAcid(letter) {
        switch (letter) {
            //Phe
            case "F": return { letter: "F", acronym: "Phe", name: "Phenylalanine"};
            //Leu
            case "L": return {letter: "L",  acronym: "Leu", name: "Leucine"};
            //Ser
            case "S": return { letter: "S", acronym: "Ser", name: "Serine"};
            //Tyr
            case "Y": return { letter: "Y", acronym: "Tyr", name: "Tyrosine"};
            //STOP
            case "*": return { letter: "*", acronym: "STOP", name: "STOP"};
            //Cys
            case "C": return { letter: "C", acronym: "Cys", name: "Cysteine"};
            //Trp
            case "W": return {letter: "W",  acronym: "Trp", name: "Tryptophan"};
            //Pro
            case "P": return { letter: "P", acronym: "Pro", name: "Proline"};
            //His
            case "H": return { letter: "H", acronym: "His", name: "Histidine"};
            //Gln
            case "Q": return {letter: "Q",  acronym: "Gln", name: "Glutamine"};
            //Arg
            case "R": return { letter: "R", acronym: "Arg", name: "Arginine"};
            //Ile
            case "I": return {letter: "I",  acronym: "Ile", name: "Isoleucine"};
            //Met
            case "M": return { letter: "M", acronym: "Met", name: "Methionine"};
            //Thr
            case "T": return {letter: "T",  acronym: "Thr", name: "Threonine"};
            //Asn
            case "N": return { letter: "N", acronym: "Asn", name: "Asparagine"};
            //Lys
            case "K": return { letter: "K", acronym: "Lys", name: "Lysine"};
            //Val
            case "V": return {letter: "V",  acronym: "Val", name: "Valine"};
            //Ala
            case "A": return {letter: "A",  acronym: "Ala", name: "Alanine"};
            //Asp
            case "D": return {letter: "D",  acronym: "Asp", name: "Aspartate"};
            //Glu
            case "E": return { letter: "E", acronym: "Glu", name: "Glutamate"};
            //Gly
            case "G": return { letter: "G", acronym: "Gly", name: "Glycine" };
            //
            default:
                return null
        }
    }

    map_AA_scores(numbers){
        const aminoAcids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*"];
        if (numbers.length != aminoAcids.length){
            return null
        }
        var iterator=0
        aminoAcids.forEach((e)=>{
            this.map.set(e,numbers[iterator]); //associate each amino acid with their blosum score
            iterator+=1;
        })
        return this.map
    }


    get_AA_score(symbol){
        return this.map.get(symbol)
    }
}



//Handles the querying of aminoacid pairs to return their score
//Built using AminoScore objects that are mapped to their symbol
class BlosumMatrix{
    constructor(){
        this.matrix=new Map(); //pairs ("aa1",AminoScore)
    }

    get_pair_score(aa1,aa2){
        return this.matrix.get(aa1).get_AA_score(aa2)
    }

    add_aminoscore(aa1,obj){
        this.matrix.set(aa1,obj)
    }

}



//Handles the reading and formatting of blosum62.txt
//
//can output the BlosumMatrix object that can be queried
class BlosumMatrixParser{

    constructor(path="./BLOSUM62.txt"){
        this.path=path
        this.data=""
        this.lines=[]
        this.scores=[]
        this.init()
    }

    init(){
    try{
        const data = fs.readFileSync(this.path, 'utf8');
        this.data = data.toString()
        this.lines = this.data.split("\n")
        //console.log("init done")
      } 
      catch (err) {
        console.error(err);
        }
    }

    readLine(i){
        return this.lines[i]
    }

    get_array_at_row_n(n){
        return  this.lines[n].split(/\s+/)
    }


    construct_matrix(){
        var blosum_matrix = new BlosumMatrix();

        for(var i =1;i<this.lines.length;i++){
            //get lines with scores & symbol
            var numbers = this.get_array_at_row_n(i)
            var symbol = numbers[0].toString()
            //construct aminoscore obj
            if(numbers[numbers.length-1]==""){
                
                numbers =  numbers.splice(1,numbers.length -2)
            }
            else{
                numbers =  numbers.splice(1,numbers.length -1)
            }
            var obj = new AminoScore(symbol)
            obj.map_AA_scores(numbers)
            //
            blosum_matrix.add_aminoscore(symbol,obj)
        }
        return blosum_matrix
    }
}



class FastaParser{
    constructor(path="./sequences.fasta"){
        this.path=path
        this.data=""
        this.lines=[]
        this.sequences=[]
        this.init()
    }

    init(){
        try{
            const data = fs.readFileSync(this.path, 'utf8');
            this.data = data.toString()
            this.lines = this.data.split("\n")
            //console.log("init done")
          } 
          catch (err) {
            console.error(err);
            }
        }
    
    readLine(i){
        return this.lines[i]
    }

    list_sequences(){
        var sequences = []
        for(var i = 0; i<this.lines.length; i++)
            if (this.lines[i].includes(">")){

                var sequence = this.lines[i+1]+this.lines[i+2]+this.lines[i+3]
                sequences.push(sequence);
                i+=3
            }
        return sequences
    }
}


/*
Rappel sur la matière sur BLOSUM

BLOCK = petit site très conservé évolutivement
BLOSUM = Procure des matrices de probablité d'échange d'un aa à un autre
BLOSUM80 = Matrice construite à partir de séquences 80% identiques


Différences par rapport à PAM:
Table générées via séquences hautement conservées basées sur le % d'identité plutôt que la fréquence de mutations / distance évolutive
Ne s'intéresse pas au lien évolutif comme PAM
On cherche un alignement local dans un sous-domaine
Les matrices sont issues des données, PAM extrapole les matrices via PAM1
BLOSUM62 est utilisé dans l'alignement BLAST


Score(aa1,aa2) = Fréquence observée de (aa1,aa2) / (Fréquence(aa1)*Fréquence(aa2))

Un score positif signifie que pour des séquences avec un minimum d'identité (%)
le changement se produit plus fréquemment que le hasard le prévoit.
*/



//-----------------------------------------Fin parsers


//Start algos------------------------------


function makeArray(rows,cols,value){
    var array = []
    for (var i = 0; i < rows; i++) {
        array.push([])
        for (var j = 0; j < cols; j++) {
            array[i].push(value);
        }
    }
    return array
}


function LevensteinDistance(seq1,seq2){
    //
    var dp = makeArray(seq1.length+1,seq2.length+1,0) //array of zeros
    for(var i = 0; i < dp.length; i++){
        dp[i][0]=i //first row = 0,1,2,3...
        if (i==0){
            for(var j=0;j<dp[i].length;j++){  //first col = 0,1,2,3...
                dp[i][j]=j
            }
        }
    }
    //dp array should be ready at this step
    //
    //recursive step :
    for (var i = 1; i < dp.length; i++) {
        for (var j = 1; j < dp[i].length; j++) {
            //
            var case1 = dp[i-1][j]+1
            var case2 = dp[i][j-1]+1
            var case3 = 0
            seq1[i]==seq2[j]?case3=dp[i-1][j-1]:case3=dp[i-1][j-1]+1

            dp[i][j] = Math.min(case1,case2,case3);
        }
    }
    return dp[dp.length-1][dp[dp.length-1].length-1] //return bottom right value
}



function distanceEdition(sequences){
    //
    var ed = makeArray(sequences.length,sequences.length,0)

    for(var i=0; i<sequences.length; i++){
        for(var j=0; j<sequences.length; j++){
            if(i==j){
                continue
            }
            ed[i][j]=LevensteinDistance(sequences[j],sequences[i])
        }
    
    }
    return ed;

}

function sequenceCentrale(sequences){
    //
    var de = distanceEdition(sequences)
    //
    var min_score = 10000000000000
    var best_seq = null
    //
    var it=0
    de.forEach(col=>{
        var score = 0
        col.forEach(value=>{score+=value})
        if (score<min_score){
            min_score=score
            best_seq=sequences[it]
        }
        it+=1
    })
    console.log(`best score ${min_score} with sequence ${best_seq} `)
    return best_seq
}

function LevensteinDistance(seq1,seq2){
    //
    var dp = makeArray(seq1.length+1,seq2.length+1,0) //array of zeros
    for(var i = 0; i < dp.length; i++){
        dp[i][0]=i //first row = 0,1,2,3...
        if (i==0){
            for(var j=0;j<dp[i].length;j++){  //first col = 0,1,2,3...
                dp[i][j]=j
            }
        }
    }
    //dp array should be ready at this step
    //
    //recursive step :
    for (var i = 1; i < dp.length; i++) {
        for (var j = 1; j < dp[i].length; j++) {
            //
            var case1 = dp[i-1][j]+1
            var case2 = dp[i][j-1]+1
            var case3 = 0
            seq1[i]==seq2[j]?case3=dp[i-1][j-1]:case3=dp[i-1][j-1]+1

            dp[i][j] = Math.min(case1,case2,case3);
        }
    }
    return dp[dp.length-1][dp[dp.length-1].length-1] //return bottom right value
}

class Alignment{



}


//Utilisation du parsers et de l'objet blosumMatrix : 
var parser = new BlosumMatrixParser("./BLOSUM62.txt")
var blosumMatrix = parser.construct_matrix()
console.log(blosumMatrix.get_pair_score("A","S"))

var parser2 = new FastaParser("./sequences.fasta")
var sequences = parser2.list_sequences()



sequenceCentrale(sequences)