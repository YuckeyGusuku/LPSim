#include "accelerationtest.h"

const int rows = 481;//行の大きさ
const int cols = 455;//列の大きさ
const int centerRow = 241;//フィルタの中心地
const int centerCol = 241;//フィルタの中心地
const int amountmode = 8; //lpモード数
double minloss = 100;
string FileName = "lp6-1O-SI-diff.txt";
int d1 = 9, d2 = 19, d3 = 51, d4 = 51, d5 = 51, d6 = 51, d7 = 51, d8 = 51;

int main() {
    const string files[amountmode] = {"lp01mode.txt", "lp11mode.txt", "lp02mode.txt", "lp21mode.txt","lp31mode.txt","lp12mode.txt"/*,"lp41mode.txt","lp22mode.txt"*/};
    vector<double> matrix(amountmode * rows * cols, 0.0);  // 光強度分布を入れる配列
    string line;
    double lppower[amountmode] = {0.0};//lpパワー,損失
    double aplppower[amountmode] = {0.0};//通過後lpパワー
    double diffloss = 0.0; //全体での最大誤差
    // ファイル読み込み
    #pragma omp parallel for
    for (int i = 0; i < amountmode; i++) {
        ifstream infile(files[i]);
        int row = 0;
        while (getline(infile, line)) {
            istringstream iss(line);
            string value;
            int col = 0;
            while (getline(iss, value, '\t')) {
                double val = stod(value);
                matrix[i * rows * cols + row * cols + col] = val;
    
                // atomic操作により並列化による競合回避
                #pragma omp atomic
                lppower[i] += val;
                col++;
            }
            row++;
        }
        infile.close();
    }

    ofstream file(FileName);
    if (!file) {
        // ファイルが開けなかった場合
        std::cerr << "Error: ファイルを開けませんでした。ファイル名: " << FileName << std::endl;
        return 1; // エラー終了
    }

    //時間測定
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    LARGE_INTEGER start, end;
    QueryPerformanceCounter(&start);
    createfilter(matrix,aplppower);
    calculateloss(lppower, aplppower, file);
    //cout << lppower[0]<<endl;
    //フィルタの作成と通過後のlpパワー測定
    // #pragma omp parallel for schedule(dynamic)
    // for(d8 = 50; d8 > 4 ; d8--){
    //    for(d7 = 49; d7 >= d6; d7--){
    //        for(d6 = 49; d6 > d5; d6--){
    //            for(d5 = 48;d5 >= d4; d5--){
    //                for(d4 = 48;d4 > d3; d4--){
    //                    for(d3 = 47; d3 >= d2; d3--){
    //                        for(d2 = 47; d2 > d1; d2--){
    //                                for(d1 = 46 ;d1 >=1 ; d1--){
    //                                    createfilter(matrix,aplppower);
    //                                    #pragma omp critical
    //                                    {
    //                                        calculateloss(lppower, aplppower, file);
    //                                    }
    //                                }
    //                        }
    //                    }
    //                }
    //            }
               
    //        }
    //    }
    //    cout <<"\n"<<"------------------------d8 = "<< d8<< "\n";
    //    file.flush();
    // }

    //フィルタの作成と通過後のlpパワー測定
    // #pragma omp parallel for schedule(dynamic)
    // for(d8 = 50; d8 > 4 ; d8--){
    //    for(d7 = 49; d7 >= d6; d7--){
    //        for(d6 = 49; d6 > d5; d6--){
    //            for(d5 = 48;d5 >= d4; d5--){
                    //for(d4 = 60;d4 > 4; d4--){
                        //for(d3 = 60; d3 > 3; d3--){
                           //for(d2 = 60; d2 > 2; d2--){
                                    for(d1 = 60 ;d1 >=1 ; d1--){
                                        createfilter(matrix,aplppower);
                                        #pragma omp critical
                                        {
                                            calculateloss(lppower, aplppower, file);
                                        }
                                   }
                           //}
                       //}
                   //}
    //            }
               
    //        }
    //    }
    //    cout <<"\n"<<"------------------------d8 = "<< d8<< "\n";
    //    file.flush();
    // }
    QueryPerformanceCounter(&end);
    double spead = static_cast<double>(end.QuadPart - start.QuadPart) * 1000.000 / freq.QuadPart;
    
    cout << "time=" << spead << "msec\n";
    
    cout << "Press Enter to continue...";
    cin.ignore();
    file.close();
    return 0;

}

double calculateValue(int row, int col) {
    double distance = sqrt(pow(row - centerRow, 2) + pow(col - centerCol, 2));
    
    if(distance <= d1){
        return 0;
    }else{
        return 1;
    }
    // Calculate value based on distance
    /*if(distance <= d3){
        return 0;
    }else if(d3 < distance && distance < d4){
        return  (distance - d3) / (d4 - d3);
    }else if (d4 <= distance && distance <= d5) {
        return 1;
    } else if (d5 < distance && distance < d6) {
        return 1-(distance - d5) / (d6 - d5);
    } else if (d6 <= distance && distance <= d7) {
        return 0;
    } else if (d7 < distance && distance < d8) {
        return  (distance - d7) / (d8 - d7);
    } else if (d6 <= distance && distance <= d7) {
        return 0;
    } else if (d7 < distance && distance < d8) {
        return (distance - d7) / (d8 - d7);
    }*else if (d8 <= distance) {
        return 1;
    } else {
        return 0;
    }*/
}

double createfilter(vector<double>& matrix , double *aplppower) {
    // フィルタの作成と計算
    const int blockSize = 100; // キャッシュブロッキング
    int j1 = 0, j2 = 0, j3 = 0, j4 = 0;
    
    
    #pragma omp parallel for collapse(3) schedule(dynamic) private(j1, j2, j3, j4)
    for (int n = 0; n < amountmode; n ++) {
        for (int bi = 0; bi < rows; bi += blockSize) {
            for (int bj = 0; bj < cols; bj += blockSize) {
                for (int i = bi; i < min(cols, bi + blockSize); i++) {
                    for (int j = bj; j < min(rows,bj + blockSize); j += 5) {
                            j1 = j + 1;
                            j2 = j + 2;
                            j3 = j + 3;
                            j4 = j + 4;
                            double gridVal = calculateValue(i, j) * matrix[n * rows * cols + i * cols + j];
                            double gridVal1 = calculateValue(i, j1) * matrix[n * rows * cols + i * cols + j1];
                            double gridVal2 = calculateValue(i, j2) * matrix[n * rows * cols + i * cols + j2];
                            double gridVal3 = calculateValue(i, j3) * matrix[n * rows * cols + i * cols + j3];
                            double gridVal4 = calculateValue(i, j4) * matrix[n * rows * cols + i * cols + j4];

                            #pragma omp atomic
                            aplppower[n] += gridVal + gridVal1 + gridVal2 + gridVal3 + gridVal4;
                    }
                }
            }
        }
    }
    //cout << "pplppwoer = "<<aplppower[0] << endl;

    return 0;
}

double calculateloss(double *lppower, double *aplppower, ofstream& file){
    double lploss[amountmode] = {0.0};
    // loss計算
    //lploss[0] = 10 * log10(lppower[0] / aplppower[0]);
    //lploss[1] = 10 * log10(lppower[1] / aplppower[1]);
    //lploss[2] = 10 * log10(lppower[2] / aplppower[2]);
    //lploss[3] = 10 * log10(lppower[3] / aplppower[3]);
    //lploss[4] = 10 * log10(lppower[4] / aplppower[4]);
    //lploss[5] = 10 * log10(lppower[5] / aplppower[5]);
    for(int i = 0; i < amountmode;i++){
        lploss[i] = 10 * log10(lppower[i] / aplppower[i]);
    }
    //cout << "\n LP12loss = "<<lploss[5] <<"dB LP12power = "<< lppower[5] << "dB apLP12power = " << aplppower[5]<<"dB" << endl;

    double difflp11 = abs(lploss[0]-lploss[1]-1);
    double difflp02 = abs(lploss[0]-lploss[2]-1.8);
    double difflp21 = abs(lploss[0]-lploss[3]-1.8);
    double difflp31 = abs(lploss[0]-lploss[4]-2.6);
    double difflp12 = abs(lploss[0]-lploss[5]-2.6);
    /*
    double difflp22 = abs(lploss[0]-lploss[6]-3.4);
    double difflp41 = abs(lploss[0]-lploss[7]-3.4);*/

    //cout <<difflp11<<","<<difflp02<<","<<difflp21<<"lploss = "<<lploss[0]<<"\n";
    //file <</*max({difflp11,difflp02,difflp21})*/difflp11<<"lploss = "<<lploss[0]<<"aplppower = "<<aplppower[0]<<"\n";
    //cout <<"lp01loss="<< lploss[0]<< "d1~8="<<d1<<","<<d2<<","<<d3<<","<<d4<<","<<d5<<","<<d6<<","<<d7<<","<<d8<< "\n";
    //cout << lppower[0] <<"aplppower = "<<aplppower[0]<<"\n";
    //cout << "difflp22   "<<difflp22 <<", difflp41   "<< difflp41<<", difflp11"<< difflp11<<", difflp02"<< difflp02<<", difflp21"<< difflp21<<", difflp31"<< difflp31<<"\n";
   
    //もし損失差誤差0.1dB以内かつ、最小値である場合のみ出力
     //if(max({difflp11,difflp02,difflp21,difflp31,difflp12}) <= 0.1){
     //    cout <<"lploss="<< lploss[0]<<"minloss"<<minloss<<"\n";
     //    if(minloss >lploss[0]){    
     //        file <<"lploss01="<< lploss[0]<<"dB,lp11="<<lploss[1]<<"dB,"<<lploss[2]<<"dB,"<<lploss[3]<<"dB"<<lploss[4]<<"dB,"<<lploss[5]<<"dB\n";
     //        file<<"d1="<<d1<<","<<d2<<","<<d3<<","<<d4/*<<","<<d5<<","<<d6<<","<<d7<<","<<d8*/<< "\n";
     //        cout <<"lploss="<< lploss[0]<<"dB,lp11="<<lploss[1]<<"dB,\n"<<lploss[2]<<"dB,"<<lploss[3]<<"dB"<<lploss[4]<<"dB,"<<lploss[5]<<"dB\n";
     //        file <<"誤差lp11="<<difflp11<< "db"<<"誤差lp02="<<difflp02<< "db"<<"誤差lp21="<<difflp21<< "db"<<"誤差lp31="<<difflp31<< "db"<<"誤差lp12="<<difflp12<< "db\n";
     //        minloss = lploss[0];
     //    }
     //}
     double maxdiff = max({difflp11,difflp02,difflp21,difflp31,difflp12/*,difflp22,difflp41*/});
     //cout << "maxdiff:   " << maxdiff<<"\n";
     if( maxdiff <= minloss){
         cout <<"lploss="<< lploss[0]<<"minloss"<<minloss<<"\n";    
             //file <<"lploss01="<< lploss[0]<<"dB,lp11="<<lploss[1]<<"dB,"<<lploss[2]<<"dB,"<<lploss[3]<<"dB"<<lploss[4]<<"dB,"<<lploss[5]<<"dB\n";
             file<<"d1="<<d1<<","<<d2<<","<<d3<<","<<d4/*<<","<<d5<<","<<d6<<","<<d7<<","<<d8*/<< "\n";
             //cout <<"lploss="<< lploss[0]<<"dB,lp11="<<lploss[1]<<"dB,\n"<<lploss[2]<<"dB,"<<lploss[3]<<"dB"<<lploss[4]<<"dB,"<<lploss[5]<<"dB\n";
             cout<<"d1="<<d1<<","<<d2<<","<<d3<<","<<d4<<"\n";
             file <<"max誤差="<<maxdiff<<"\n";
             cout <<"max誤差="<<maxdiff<<"\n";
             minloss = maxdiff;
     }

    aplppower[0] = 0.00;
    aplppower[1] = 0.00;
    aplppower[2] = 0.00;
    aplppower[3] = 0.00;
    aplppower[4] = 0.00;
    aplppower[5] = 0.00;
    /*
    aplppower[6] = 0.00;
    aplppower[7] = 0.00;*/

    return 0;
}