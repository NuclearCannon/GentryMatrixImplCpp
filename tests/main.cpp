#include "flints.hpp"
#include "ntt.hpp"
#include "twisted_ntter.hpp"
#include "ziq_array.hpp"
#include <iostream>


int test_ntt_standard()
{
    printf("test_ntt_standard\n");
    // n=8 q=1601
    // root-n: 408
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4", "5", "6", "7", "8"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, root;
    fmpz_init(q);
    fmpz_init(root);

    string_to_fmpz("1601", q);
    string_to_fmpz("408", root);

    fmpz_vector b(8);

    ntt_standard_flint(a, b, root, 8, q);

    std::vector<std::string> result = b.export_to_vec_str();

    // 预期输出：[36, 1365, 156, 1045, 1597, 548, 1437, 228]
    std::vector<std::string> expected = {
        "36", "1365", "156", "1045", "1597", "548", "1437", "228"
    };

    int first_error = -1;
    for(int i=0;i<8;i++)
    {
        if(result[i] != expected[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result[i]<<std::endl;
            first_error = i;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(root);

    if (first_error == -1)
    {
        printf("test_ntt_standard pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_standard error...\n");
        return 0;
    }
}


int test_ntt_standard_with_ctx()
{
    printf("test_ntt_standard_with_ctx\n");
    // n=8 q=1601
    // root-n: 408
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4", "5", "6", "7", "8"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, root;
    fmpz_mod_ctx_t ctx;

    fmpz_init(q);
    fmpz_init(root);

    string_to_fmpz("1601", q);
    fmpz_mod_ctx_init(ctx, q);
    string_to_fmpz("408", root);

    fmpz_vector b(8);

    ntt_standard_flint(a, b, root, 8, ctx);

    std::vector<std::string> result = b.export_to_vec_str();

    // 预期输出：[36, 1365, 156, 1045, 1597, 548, 1437, 228]
    std::vector<std::string> expected = {
        "36", "1365", "156", "1045", "1597", "548", "1437", "228"
    };

    int first_error = -1;
    for(int i=0;i<8;i++)
    {
        if(result[i] != expected[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result[i]<<std::endl;
            first_error = i;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(root);

    if (first_error == -1)
    {
        printf("test_ntt_standard_with_ctx pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_standard_with_ctx error...\n");
        return 0;
    }
}



int test_ntt_xy()
{
    printf("test_ntt_xy\n");
    // n=8 q=1601
    // root-4n: 356
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4", "5", "6", "7", "8"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, zeta;
    fmpz_init(q);
    fmpz_init(zeta);

    string_to_fmpz("1601", q);
    string_to_fmpz("356", zeta);

    TwistedNtterXY ntter(
        8, q, zeta
    );

    fmpz_vector b(8), c(8);

    ntter.ntt(a, b);    // b=NTT(a)
    ntter.intt(b, c);   // c=iNTT(b)


    std::vector<std::string> result_b = b.export_to_vec_str();
    std::vector<std::string> result_c = c.export_to_vec_str();

    // 预期输出：[1100, 1246, 740, 1028, 1189, 1322, 332, 1056]
    std::vector<std::string> expected = {
        "1100", "1246", "740", "1028", "1189", "1322", "332", "1056"
    };

    bool error = 0;
    for(int i=0;i<8;i++)
    {
        if(result_b[i] != expected[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result_b[i]<<std::endl;
            error = 1;
            break;
        }
    }
    // 预期c==a
    for(int i=0;i<8;i++)
    {
        if(result_c[i] != arr_str[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<arr_str[i]<<" but got "<<result_c[i]<<std::endl;
            error = 1;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(zeta);

    if (!error)
    {
        printf("test_ntt_xy pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_xy error...\n");
        return 0;
    }
}


int test_ntt_w()
{
    printf("test_ntt_w\n");
    // p=5 q=1601 g=2
    // eta: 442
    // [20 741 683 152]
    std::vector<std::string> arr_str = {
        "1", "2", "3", "4"
    };
    fmpz_vector a(arr_str);
    fmpz_t q, eta;
    fmpz_init(q);
    fmpz_init(eta);

    string_to_fmpz("1601", q);
    string_to_fmpz("442", eta);

    TwistedNtterW ntter(
        5, 2, q, eta
    );

    fmpz_vector b(4), c(4);
    ntter.ntt(a, b);    // b=NTT(a)
    ntter.intt(b, c);   // c=iNTT(b)
    std::vector<std::string> result_b = b.export_to_vec_str();
    std::vector<std::string> result_c = c.export_to_vec_str();

    // 预期输出：[20 741 683 152]
    std::vector<std::string> expected = {"20", "741", "683", "152"};

    bool error = 0;
    for(int i=0;i<4;i++)
    {
        if(result_b[i] != expected[i])
        {
            std::cout << "1.error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result_b[i]<<std::endl;
            error = 1;
            break;
        }
    }
    // 预期c==a
    for(int i=0;i<4;i++)
    {
        if(result_c[i] != arr_str[i])
        {
            std::cout << "2.error at idx="<< i << ", expect "<<arr_str[i]<<" but got "<<result_c[i]<<std::endl;
            error = 1;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(eta);

    if (!error)
    {
        printf("test_ntt_w pass!\n");
        return 1;
    }
    else
    {
        printf("test_ntt_w error...\n");
        return 0;
    }
}

int test_ziq_ctx()
{
    printf("test_ziq_ctx\n");
    // n=8 p=5 q=1601
    // zeta: 356  eta: 442
    int n=8, p=5, g=2;
    int size_ = 2*(p-1)*n*n;
    fmpz_t q, zeta, eta;
    fmpz_init(q);
    fmpz_init(zeta);
    fmpz_init(eta);

    string_to_fmpz("1601", q);
    string_to_fmpz("356", zeta);
    string_to_fmpz("442", eta);
    ZiqArrayContext ctx(n, p, g, q, zeta, eta);
    fmpz_vector data(size_);
    for(int i=0;i<size_;i++)fmpz_set_ui(data[i], i);

    int error0 = -1, error1 = -1;

    auto r1 = ctx.iw_ntt(data);
    auto r2 = ctx.iw_intt(r1);
    for(int i=0;i<size_;i++)
    {
        if(fmpz_cmp(r2[i], data[i]) != 0)
        {
            error0 = i;
            break;
        }
    }

    auto r3 = ctx.xy_ntt(data);
    auto r4 = ctx.xy_intt(r3);
    for(int i=0;i<size_;i++)
    {
        if(fmpz_cmp(r4[i], data[i]) != 0)
        {
            error1 = i;
            break;
        }
    }


    fmpz_clear(q);
    fmpz_clear(zeta);
    fmpz_clear(eta);

    if ((error0 == -1) && (error1 == -1))
    {
        printf("test_ziq_ctx pass!\n");
        return 1;
    }
    else
    {
        printf("test_ziq_ctx error(%d, %d)\n", error0, error1);
        return 0;
    }
}

int test_circledast()
{
    printf("test_circledast\n");
    // n=8 p=5 q=1601
    // zeta: 356  eta: 442
    int n=8, p=5, g=2;
    int size_ = 2*(p-1)*n*n;
    fmpz_t q, zeta, eta;
    fmpz_init(q);
    fmpz_init(zeta);
    fmpz_init(eta);

    string_to_fmpz("1601", q);
    string_to_fmpz("356", zeta);
    string_to_fmpz("442", eta);
    ZiqArrayContext ctx(n, p, g, q, zeta, eta);
    fmpz_vector data(size_);
    for(int i=0;i<size_;i++)fmpz_set_ui(data[i], i);

    fmpz_vector result = ctx.circledast(data, data);
    auto result_str = result.export_to_vec_str();

    // 

    std::vector<std::string> expected = {
        "1477", "100", "324", "548", "772", "996", "1220", "1444", "1555", "690", 
        "1426", "561", "1297", "432", "1168", "303", "32", "1280", "927", "574", 
        "221", "1469", "1116", "763", "110", "269", "428", "587", "746", "905", 
        "1064", "1223", "188", "859", "1530", "600", "1271", "341", "1012", "82", 
        "266", "1449", "1031", "613", "195", "1378", "960", "542", "344", "438", 
        "532", "626", "720", "814", "908", "1002", "422", "1028", "33", "639", 
        "1245", "250", "856", "1462", "1162", "679", "196", "1314", "831", "348", 
        "1466", "983", "346", "375", "404", "433", "462", "491", "520", "549", 
        "1131", "71", "612", "1153", "93", "634", "1175", "115", "315", "1368", 
        "820", "272", "1325", "777", "229", "1282", "1100", "1064", "1028", "992", 
        "956", "920", "884", "848", "284", "760", "1236", "111", "587", "1063", 
        "1539", "414", "1069", "456", "1444", "831", "218", "1206", "593", "1581", 
        "253", "152", "51", "1551", "1450", "1349", "1248", "1147", "952", "1363", 
        "173", "584", "995", "1406", "216", "627", "843", "165", "1088", "410", "1333", 
        "655", "1578", "900", "734", "568", "402", "236", "70", "1505", "1339", "1173", 
        "625", "971", "1317", "62", "408", "754", "1100", "1446", "516", "1374", "631", 
        "1489", "746", "3", "861", "118", "407", "176", "1546", "1315", "1084", "853", 
        "622", "391", "298", "579", "860", "1141", "1422", "102", "383", "664", "189", 
        "982", "174", "967", "159", "952", "144", "937", "847", "551", "255", "1560", 
        "1264", "968", "672", "376", "1445", "60", "276", "492", "708", "924", "1140", 
        "1356", "442", "1170", "297", "1025", "152", "880", "7", "735", "1040", "679", 
        "318", "1558", "1197", "836", "475", "114", "37", "188", "339", "490", "641", 
        "792", "943", "1094", "635", "1298", "360", "1023", "85", "748", "1411", "473",
        "1233", "807", "381", "1556", "1130", "704", "278", "1453", "230", "316", "402", 
        "488", "574", "660", "746", "832", "847", "1445", "442", "1040", "37", "635", 
        "1233", "230", "551", "60", "1170", "679", "188", "1298", "807", "316", "255", 
        "276", "297", "318", "339", "360", "381", "402", "1560", "492", "1025", "1558", 
        "490", "1023", "1556", "488", "1264", "708", "152", "1197", "641", "85", "1130", 
        "574", "968", "924", "880", "836", "792", "748", "704", "660", "672", "1140", "7", 
        "475", "943", "1411", "278", "746", "376", "1356", "735", "114", "1094", "473", 
        "1453", "832", "952", "843", "734", "625", "516", "407", "298", "189", "1363", 
        "165", "568", "971", "1374", "176", "579", "982", "173", "1088", "402", "1317", 
        "631", "1546", "860", "174", "584", "410", "236", "62", "1489", "1315", "1141", 
        "967", "995", "1333", "70", "408", "746", "1084", "1422", "159", "1406", "655", 
        "1505", "754", "3", "853", "102", "952", "216", "1578", "1339", "1100", "861", 
        "622", "383", "144", "627", "900", "1173", "1446", "118", "391", "664", "937", 
        "1162", "346", "1131", "315", "1100", "284", "1069", "253", "679", "375", "71", 
        "1368", "1064", "760", "456", "152", "196", "404", "612", "820", "1028", "1236", 
        "1444", "51", "1314", "433", "1153", "272", "992", "111", "831", "1551", "831", 
        "462", "93", "1325", "956", "587", "218", "1450", "348", "491", "634", "777", 
        "920", "1063", "1206", "1349", "1466", "520", "1175", "229", "884", "1539", "593", 
        "1248", "983", "549", "115", "1282", "848", "414", "1581", "1147", "1477", "1555", 
        "32", "110", "188", "266", "344", "422", "100", "690", "1280", "269", "859", "1449", 
        "438", "1028", "324", "1426", "927", "428", "1530", "1031", "532", "33", "548", "561", 
        "574", "587", "600", "613", "626", "639", "772", "1297", "221", "746", "1271", "195", 
        "720", "1245", "996", "432", "1469", "905", "341", "1378", "814", "250", "1220", 
        "1168", "1116", "1064", "1012", "960", "908", "856", "1444", "303", "763", "1223", 
        "82", "542", "1002", "1462"
    };


    int error = -1;
    for(int i=0;i<size_;i++)
    {
        if (result_str[i] != expected[i])
        {
            std::cout << "error at idx="<< i << ", expect "<<expected[i]<<" but got "<<result_str[i]<<std::endl;
            error = i;
            break;
        }
    }

    fmpz_clear(q);
    fmpz_clear(zeta);
    fmpz_clear(eta);

    if (error == -1)
    {
        printf("test_circledast pass!\n");
        return 1;
    }
    else
    {
        printf("test_circledast error...\n");
        return 0;
    }
}


int main()
{
    int sum = 0;
    sum += test_ntt_standard();
    sum += test_ntt_standard_with_ctx();
    sum += test_ntt_xy();
    sum += test_ntt_w();
    sum += test_ziq_ctx();
    sum += test_circledast();

    std::cout << "Total Pass: " << sum << std::endl; 
}