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
            //console.log(numbers)
            //console.log(aminoAcids)
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
            numbers =  numbers.splice(1,numbers.length -2)
            var obj = new AminoScore(symbol)
            obj.map_AA_scores(numbers)
            //
            blosum_matrix.add_aminoscore(symbol,obj)
        }
        return blosum_matrix
    }
}


//Utilisation du parsers et de l'objet blosumMatrix : 
var parser = new BlosumMatrixParser("./BLOSUM62.txt")
var blosumMatrix = parser.construct_matrix()
console.log(blosumMatrix.get_pair_score("V","S"))