# CG_SHOP_2025_1

1. Τίτλος Project:
Μη Αμβλυγώνια Τριγωνοποίηση Επίπεδων Γράφων Ευθύγραμμων Τμημάτων (Planar Straight Line Graphs) με τη χρήση της βιβλιοθήκης CGAL.

2. Περιγραφή:
    Η μείωση των αμβλυγώνιων τριγώνων γίνεται με τις εξής συναρτήσεις:
    void start_the_flips(Custom_CDT& cdt, const Polygon& polygon);
    void insert_circumcenter_centroid(Custom_CDT& custom_cdt, const Polygon& polygon);
    void insert_projection(Custom_CDT& custom_cdt, const Polygon polygon);
    void insert_midpoint(Custom_CDT& custom_cdt, const Polygon& polygon);
    void insert_orthocenter(Custom_CDT& custom_cdt, const Polygon& polygon);

    Γενικώς οι υλοποιήσεις των παραπάνω συναρτήσεων έχουν την φιλοσοφία να προσομιάζουμε (simulate) μια εισαγωγή steiner point,
    έχοντας ένα αντίγραφο της βασικής τριγωνοποίσης (simulate cdt) και εισάγοντας steiner points σε ένα simulate cdt βλέπουμε αν
    έχει όφελος (μείωση αμβλειών) η εισαγωγή αυτή, ώστε να την εκτελέσουμε και στην βασική τριγωνοποίηση (cdt).

    Η συνάρτηση start_the_flips() διατρέχει όλες τις ακμές του cdt βρίσκοντας τα 2 faces που χωρίζονται από αυτή της ακμή
    και από αυτά τα 2 faces παίρνει τα 4 σημεία τους (Points) και ελέγχει αρχικά αν τα 4 σημεία σχηματίζουν κυρτό πολύγωνο
    μέσω της συνάρτησης can_flip() ή οποία καλεί την is_convex() για να ελέγει την κυρτότητα του πολυγώνου.
    Η συνάρτηση can_flip() στο τέλος θα μας επιστρέψει γενικά αν μας συμφέρει (λιγότερες αμβλείες) να κάνουμε flip μια ακμή
    και όχι αν απλά γίνεται flip η ακμή.

    Η συνάρτηση insert_circumcenter_centroid() εισάγει circumcenter point σε μια τριγωνοποίηση αν το circumcenter point είναι
    εντός boundary, και το τετράπλευρο που σχηματίζουν τα σημεία του face που διατρέχουμε με το circumcenter point αν είναι κυρτό,
    ώστε μετά την εισαγωγή του circumcenter point να ελέγχουμε με flips αν υπάρχει μείωση αμβλειών ώστε να το προσθέσουμε στην τριγωνοποίηση.
    Αν το circumcenter point πέφτει εκτός boundary, τότε καλείτε η can_insert_centroid() για να ελέγξει αν μας συμφέρει να εισάγουμε
    centroid point, όπου και εκεί μετά την εισαγωγή centroid χρησιμοποιούμε την start_the_flips() για να δούμε αν και με centroid point και με flips μειώνουμε τις αμβλείες της τριγωνοποίσης.

    Η συνάρτηση insert_projection() διατρέχει τα faces του cdt, όταν βρεί αμβλυγώνιο τρίγωνο, βρίσκει σε ποιο point του face συμβαίνει αυτό
    αλλά και την ακμή απέναντι από την αμβλεία γωνία, ώστε να καλέσει την συνάρτηση line.projection() για να πάρουμε το σημείο προβολής της 
    αμβλείας γωνίας με την απέναντι ευθεία.
    (1) Σαν πρώτο βήμα διατρέχουμε όλα τα faces της cdt αλλά ελέγχουμε μόνο τα τρίγωνα που η μια τους ακμή είναι πάνω στο region boundary και εισάγουμε το projection point αν δούμε μείωση αμβλειών.
    Ύστερα ακολουθούν δύο περιπτώσεις κατά την νέα διάσχιση των faces του cdt.
    Η πρώτη απλώς ελέγχει αν η εισαγωγή projection point μειώνε τα αμβυγώνια και αν ναι προσθέτει το projection point στο cdt.
    Η δεύτερη (και πιο σημαντική) περίπτωση είναι αν δεν υπάρχει βελτίωση στην μείωση των αμβλυγωνίων αλλά ούτε και αύξηση μετά την εισαγωγή projection point και αν το projection point (κατά την προσομείωση του insert) πλέον είναι point ενός τριγώνου που αυτό το τρίγωνο έχει ακμή στο region boundary και το projection point είναι το εσωτερικό point αυτού του τριγώνου (δεν είναι κάποιο από τα points της ακμής που είναι πάνω στο region boundary), τότε αυτό το point το εισάγουμε (αν και προσωρινά δεν έχουμε όφελος) και καλούμε την περίπτωση (1) που θα διατρέξει και θα ελέγξει μόνο τα faces που η ακμή τους είναι πάνω στο region boundary, ώστε μετά απλώς να προστεθεί μια προβολή από το projection point στην ακμή του region boundary και έτσι τα αμβλυγώνια θα μειωθούν κατά 1.

    Η συνάρτηση insert_midpoint() και η συνάρτηση insert_orthocenter() εισάγουν στο cdt ένα midpoint και ένα orthocenter point εφόσον
    πάλι η εισαγωγή τους μας μειώνει τα αμβλυγώνια, κάνοντας παράλληλα και τα κατάλληλα flips μετά την εισαγωγή τους.

    Στο τέλος μετρούνται τα steiner points που έγιναν εισαγωγή και οι ακμές που δημιούργησαν και αποθηκεύονται στο solution_output.json
    που δημιουργείτε στο τέλος του προγράμματος.

    Δεν χρησιμοποιήθηκε εξαντλητική μέθοδος των παραπάνω συναρτήσεων αλλά η κλήση τους γίνεται ατομικά για να μειώσουμε τον χρόνο επιστροφής
    των αποτελεσμάτων. Αυτό σημαίνει ότι για κάθε cdt δεν είναι τα αποτελέσματα βέλτιστα, αλλά σίγουρα χρησιμοποιήθηκαν λιγότερα steiner points
    και η επιστροφή των αποτελεσμάτων είναι ακριαία.

    Υπάχρουν και όλα τα tests που δόθηκαν, στον φάκελο tests ώστε να ελεγθεί το πρόγραμμα, και στο tests/test_doc.txt βρίσκονται έτοιμες 
    οι εντολές για να διαβάσετε ένα json test και να προσθέσετε απλά την εντολή στην γραμμή 17 του project.cpp.

    Το link του Github Repository είναι το εξής:
    https://github.com/MikelosR/CG_SHOP_2025_1

3. Οργάνωση Φακέλων:
CG_SHOP_2025_1:
    /tests: 
        a. Τα .json instances που δόθηκαν στο eclass ώστε να ελέγξουμε τον κώδικά μας.
        b. Aρχείο test_doc.txt που καταγράφουμε για κάθε τεστ που δόθηκε το ποσοστό επιτυχίας της τριγωνοποίησης.

    /includes/utils: 
        a. Custom_Constrained_Delaunay_triangulation_2.h όπως δόθηκε στις συζητήσεις του eclass του μαθήματος.
        b. functions.cpp οι υλοποιήσεις των συναρτήσεων που δημιουργήσαμε ώστε να εκπονήσουμε την εργασία.
        c. functions.h οι δηλώσεις των συναρτήσεων
    
    CMakeLists.txt: 
        Το αρχείο που χρειάζεται για την εισαγωγή των κατάλληλων βιβλιοθηκών (CGAL, Boost, γραφικών) ώστε να παράγουμε το επιθυμητό εκτελέσιμο.
    
    project.cpp:
        Το αρχείο μας με την main function που αντλεί δεδομένα από ένα .json αρχείο με δεδομένα για έναν γράφο πάνω στον οποίο δημιουργούμε την  τριγωνοποίηση Delaunay, την βελτιστοποιούμε αποσκοπώντας στο να ελαχιστοποιήσουμε τα αμβλυγώνια τρίγωνα, με τις εξής μεθόδους:
            a. flips σε ακμές τριγώνων.
            b. εισαγωγή steiner points σε: circumcenter/centroid, midpoint, projection & orthocenter.
        Τρέχουμε επαναληπτικά τις μεθόδους flip, circumcenter - centroid, midpoint, projection, orhtocenter μέχρι να σταμτήσουν να έχουν 
        αποτέλεσμα (μείωση αμβλυών) ή μέχρι να καταφέρουν να μηδενίσουν τις αμβλείες του cdt.
        Τέλος, δημιουργούμε ένα αρχείο solution_output.json με τα αποτελέσματα όπως ζητήθηκαν από την εκφώνηση. 

    Makefile:
    Το αρχείο που δημιουργήθηκε όταν εκτελέσαμε στο terminal την εντολή: 
    cmake -DCGAL_DIR=/usr/lib/CGAL
    
    solution_output.json: 
        Το παραγόμενο output από την main project.cpp

    .git
        Ο φάκελος με τις πληροφορίες για push, commite etc στο Github Repository

    ΣΗΜΕΙΩΣΗ: Οι παρακάτω εντολές πρέπει να βρίσκονται στο αρχείο .bashrc:
        export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
        export DISPLAY=:0
        export BOOST_ROOT= "το path που βρίσκεται ο φάκελος boost_1_80_0"
        export CGAL_ROOT="το path που βρίσκεται ο φάκελος της βιβλιοθήκης CGAL-5.6.1"
    Χρησιμοποιήθηκε η έκδοση 1_80_0 της Boost στην υπολοποίησή μας.

4. Οδηγίες Μεταγλώττισης & Εκτέλεσης στο terminal:
    Step 1: /path/to/program
    Step 2: cmake -DCGAL_DIR=/usr/lib/CGAL
    Step 3: make
    Step 4: ./project

    ΠΡΟΣΟΧΗ! Δεν πρέπει να τρέξει η εντολή cgal_create_CMakeLists -s, το CMakeLists.txt στον υπάρχων φάκελο είναι το κατάλληλο
    και περιέχει τα κατάλληλα includes βιβλιοθηκών για το υπάρχον πρόγραμμα

    Αποτέλεσμα: 
        α. Εκτυπώνει κάποια μηνύματα που περιγράφουν πως επηρέασε τα αμβλυγώνια τρίγωνα η κάθε τεχνική που χρησιμοποιήσαμε.
        β. Δημιουργεί το αρχείο solution_output.json με τα αποτελέσματα όπως ζητήθηκαν στην εκφώνηση της εργασίας.

    Optional Step: 
        Για εκτέλεση πάνω σε διαφορετικά inputs, διαλέγω μία εντολή από τον φάκελο tests/test_doc.txt 
        και την αντικαθιστώ με την εντολή στη γραμμή 17 του project.cpp -> make clean -> Step 2 -> Step 3.

5. Στοιχεία Φοιτητών:
ΑΝΑΣΤΑΣΟΠΟΥΛΟΣ ΑΝΔΡΕΑΣ sdi1900009
ΠΑΠΑΔΗΜΟΠΟΥΛΟΣ ΜΙΧΑΗΛ-ΑΓΓΕΛΟΣ sdi2000163