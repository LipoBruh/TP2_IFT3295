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



class Alignement{
    constructor(score,seq1,seq2){
        this.score = score
        this.original_seq1=seq1
        this.original_seq2=seq2
        this.seq1=""
        this.seq2=""

    }


    traceback(M, X, Y, gapOpen, gapExtend, blosumMatrix) {
        let i = this.original_seq1.length;
        let j = this.original_seq2.length;
    
        // Determine the starting matrix based on the maximum score at (n, m)
        let currentMatrix = "M";
        let maxScore = Math.max(M[i][j], X[i][j], Y[i][j]);
        if (maxScore === X[i][j]) currentMatrix = "X";
        if (maxScore === Y[i][j]) currentMatrix = "Y";
    
        while (i > 0 || j > 0) {
            if (currentMatrix === "M") {
                // match/Mismatch
                const matchScore = Number(blosumMatrix.get_pair_score(this.original_seq1[i - 1], this.original_seq2[j - 1]));
    
                if (i > 0 && j > 0 && M[i][j] === M[i - 1][j - 1] + matchScore) {
                    // No gap
                    this.seq1 = this.original_seq1[i - 1] + this.seq1;
                    this.seq2 = this.original_seq2[j - 1] + this.seq2;
                    i--;
                    j--;
                } else if (i > 0 && j > 0 && M[i][j] === X[i - 1][j - 1] + matchScore) {
                    currentMatrix = "X";
                } else if (i > 0 && j > 0 && M[i][j] === Y[i - 1][j - 1] + matchScore) {
                    currentMatrix = "Y";
                }
            } else if (currentMatrix === "X") {
                // gap in seq2
                if (j > 0 && X[i][j] === M[i][j - 1] - gapOpen) {
                    currentMatrix = "M";
                } else if (j > 0 && X[i][j] === X[i][j - 1] - gapExtend) {
                }
                this.seq1 = "-" + this.seq1;
                this.seq2 = this.original_seq2[j - 1] + this.seq2;
                j--;
            } else if (currentMatrix === "Y") {
                // gap in seq1
                if (i > 0 && Y[i][j] === M[i - 1][j] - gapOpen) {
                    currentMatrix = "M";
                } else if (i > 0 && Y[i][j] === Y[i - 1][j] - gapExtend) {
                }
                this.seq1 = this.original_seq1[i - 1] + this.seq1;
                this.seq2 = "-" + this.seq2;
                i--;
            }
    
            // Handle edge cases for the first row or column
            if (i === 0 && j > 0) {
                this.seq1 = "-".repeat(j) + this.seq1;
                this.seq2 = this.original_seq2.slice(0, j) + this.seq2;
                break;
            }
            if (j === 0 && i > 0) {
                this.seq1 = this.original_seq1.slice(0, i) + this.seq1;
                this.seq2 = "-".repeat(i) + this.seq2;
                break;
            }
        }
    
        // Log the final alignment
        this.log_alignment();
    }
    


    log_alignment(){
        const reference = `${this.original_seq1}\n${this.original_seq2}\n\n\n`
        const alignment = `${this.seq1}\n${this.seq2}\nscore: ${this.score}`;
        const content = reference+alignment
        // Write to a file
        fs.writeFile('output.txt', content, (err) => {
            if (err) {
                console.error('Error writing to file:', err);
            }
        });

    }
}

class MultipleAlignment{
    constructor(){
        this.aligned_sequences=[]
    }
}




// this video https://www.youtube.com/watch?v=DQQ_q2dn2ds gives a great example
// on how to incoporate affine gap penalties to needlemanWunch
function needlemanWunschAffine(s1, s2, blosumMatrix, gapOpen=10, gapExtend=1) {
    const n = s1.length;
    const m = s2.length;

    // Initialize matrices
    const M = makeArray(s1.length+1,s2.length+1,-Infinity) 
    const X = makeArray(s1.length+1,s2.length+1,-Infinity)
    const Y = makeArray(s1.length+1,s2.length+1,-Infinity)


    // Initialize base cases
    M[0][0] = 0;
    for (let i = 1; i <= n; i++) {
        X[i][0] = -gapOpen - (i - 1) * gapExtend;
    }
    for (let j = 1; j <= m; j++) {
        Y[0][j] = -gapOpen - (j - 1) * gapExtend;
    }


    // Fill using the reccurence formula from https://www.youtube.com/watch?v=DQQ_q2dn2ds
    for (let i = 1; i <= n; i++) {
        for (let j = 1; j <= m; j++) {
            const matchScore = Number(blosumMatrix.get_pair_score(s1[i - 1], s2[j - 1]))

            // Calculate M(i, j)
            M[i][j] = Math.max(
                M[i - 1][j - 1] + matchScore,
                X[i - 1][j - 1] + matchScore,
                Y[i - 1][j - 1] + matchScore
            );
            // Calculate X(i, j)
            X[i][j] = Math.max(
                M[i - 1][j] - gapOpen,
                X[i - 1][j] - gapExtend
            );

            // Calculate Y(i, j)
            Y[i][j] = Math.max(
                M[i][j - 1] - gapOpen,
                Y[i][j - 1] - gapExtend
            );
        }
        
    }

    var alignement = new Alignement(Math.max(M[n][m], X[n][m], Y[n][m]),s1,s2)
    alignement.traceback(M,X,Y,gapOpen,gapExtend,blosumMatrix)
    return alignement
}



//Utilisation du parsers et de l'objet blosumMatrix : 
var parser = new BlosumMatrixParser("./BLOSUM62.txt")
var blosumMatrix = parser.construct_matrix()
console.log(blosumMatrix.get_pair_score("A","S"))

var parser2 = new FastaParser("./sequences.fasta")
var sequences = parser2.list_sequences()



sequenceCentrale(sequences)

var alignment = needlemanWunschAffine(sequences[0],sequences[1],blosumMatrix,10,1)
console.log("...")
console.log(alignment.seq1)
console.log(alignment.seq2)