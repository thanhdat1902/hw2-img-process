#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QPixmap>
#include <QImage>
#include <QVector>
#include <QQueue>
#include <iostream>
#include <utility>

using namespace std;
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->lineEdit_3->setText("3");
    ui->lineEdit_4->setText("3");
    ui->lineEdit_5->setText("1");
    ui->lineEdit_6->setText("3");
    ui->lineEdit_7->setText("3");

    ui->lineEdit_Q->setText("1");

    ui->lineEdit_d->setText("2");

}

MainWindow::~MainWindow()
{
    delete ui;
}
QImage imgGray;
QVector<QVector<int>> imgArray;


void MainWindow::on_pushButton_clicked()
{
    QString file_name = QFileDialog::getOpenFileName(this, tr("Open File"), QDir::homePath(), tr("Images (*.png *.xpm *.jpg *.ppm *.tif)"));
    if (!file_name.isEmpty()){

        //open prompt and display image
        QMessageBox::information(this, "...", file_name);
        QImage img(file_name);
        imgGray = img.convertToFormat(QImage::Format_Grayscale8);
        QPixmap pix = QPixmap::fromImage(imgGray);
        int w = imgGray.size().width();
        int h = imgGray.size().height();
        ui->label_pic->setSizePolicy( QSizePolicy::Ignored, QSizePolicy::Ignored );
        ui->label_pic->setPixmap(pix);

        ui->lineEdit->setText(QString::number(w));
        ui->lineEdit_2->setText(QString::number(h));

        for (unsigned int i = 0; i < w; i++){
            imgArray.push_back(QVector<int>());
            for (unsigned int j = 0; j < h; j++){
                QColor clrCurrent( imgGray.pixel( i, j ));
                int r = (clrCurrent.red() + clrCurrent.green() + clrCurrent.blue()) / 3;
                imgArray[i].push_back(r);
            }
            cout << endl;
        }
    }

}

void grayscaleToQPixmap(const QVector<QVector<int>> &grayscaleImage, QPixmap &qPixmap) {
    // Convert the grayscale image to a QImage.
    QImage image(grayscaleImage.size(), grayscaleImage[0].size(),QImage::Format_Grayscale8);
    for (int i = 0; i < grayscaleImage.size(); ++i) {
        for (int j = 0; j < grayscaleImage[i].size(); ++j) {
            QRgb value = qRgb(grayscaleImage[i][j],grayscaleImage[i][j],grayscaleImage[i][j]);
            image.setPixel(i,j,value);
        }
    }

    // Convert the QImage to a QPixmap.
    qPixmap = QPixmap::fromImage(image);
}
void MainWindow::handleRenderImage(QVector<QVector<int>> newImg) {
    // Showing result image
    int destWidth = ui->lineEdit->text().toInt();
    int destHeight = ui->lineEdit_2->text().toInt();
    QPixmap result;
    grayscaleToQPixmap(newImg , result);
    ui->label_pic_2->setPixmap(result);
    if(destWidth < 100 && destHeight < 100) {
        ui->label_pic_3->setPixmap(result.scaled(ui->label_pic_3->width(),ui->label_pic_3->height()));
    }

    ui->label_pic->adjustSize();
    // Showing new size
    ui->label_new_w->setText(QString::number(newImg.size()));
    ui->label_new_h->setText(QString::number(newImg[0].size()));

    // Reset stage
    for(int i=0 ; i< newImg.size(); i++) {
        for(int j=0 ; j< newImg[i].size(); j++) {
            newImg[i].pop_back();
        }
        newImg.pop_front();
    }
    //    imgArray = newImg;

    //    QString fileSave = "/Users/datnguyen/Downloads/ImgResult/Task2/result" + QString::number(rand() % 1000)+".png";
    //    result.save(fileSave, "PNG");
}

void getPixel(QVector<QVector<int>> img, int x, int y, unsigned char *R)
{
    // Get the colour at pixel x,y in the image and return it using the provided RGB pointers
    // Requires the image size along the x direction!
    *(R)=img[x][y];
}

double step_x,step_y;      // Step increase as per instructions above
unsigned char R1,R2,R3,R4; // Colours at the four neighbours
double RT1;                // Interpolated colours at T1 and T2
double RT2;
unsigned char R;           // Final colour at a destination pixel
int x,y;                   // Coordinates on destination image
double fx,fy;              // Corresponding coordinates on source image
double dx,dy;              // Fractional component of source image    coordinates

void handleLinearInterpolationX(QVector<QVector<int>> src,QVector<QVector<int>> &dst, int src_x, int src_y, int dest_x, int dest_y) {
    step_x=(double)(src_x-1)/(double)(dest_x-1);
    step_y=(double)(src_y-1)/(double)(dest_y-1);

    for (x=0;x<dest_x;x++) {        // Loop over destination image
        dst.push_back(QVector<int>());
        for (y=0;y<dest_y;y++)
        {
            fx=x*step_x;
            fy=y*step_y;
            dx=fx-(int)fx;
            dy=fy-(int)fy;
            getPixel(src,floor(fx),floor(fy),&R1);      // get N1 colours
            getPixel(src,ceil(fx),floor(fy),&R2);       // get N2 colours
            // Interpolate to get T1 and T2 colours
            RT1=(dx*R2)+(1-dx)*R1;

            // Store the final colour
            dst[x].push_back(RT1);
        }
    }
}

void handleLinearInterpolationY(QVector<QVector<int>> src,QVector<QVector<int>> &dst, int src_x, int src_y, int dest_x, int dest_y) {
    step_x=(double)(src_x-1)/(double)(dest_x-1);
    step_y=(double)(src_y-1)/(double)(dest_y-1);

    for (x=0;x<dest_x;x++) {        // Loop over destination image
        dst.push_back(QVector<int>());
        for (y=0;y<dest_y;y++)
        {
            fx=x*step_x;
            fy=y*step_y;
            dx=fx-(int)fx;
            dy=fy-(int)fy;
            getPixel(src,floor(fx),floor(fy),&R1);    // get N1 colours
            getPixel(src,floor(fx),ceil(fy),&R3); // get N3 colours
            // Interpolate to get T1 and T2 colours
            RT1=(dy*R3)+(1-dy)*R1;

            // Store the final colour
            dst[x].push_back(RT1);
        }
    }
}


void handleBilinearInterpolation( QVector<QVector<int>> src,QVector<QVector<int>> &dst, int src_x, int src_y, int dest_x, int dest_y)
{

    step_x=(double)(src_x-1)/(double)(dest_x-1);
    step_y=(double)(src_y-1)/(double)(dest_y-1);

    for (x=0;x<dest_x;x++) {        // Loop over destination image
        dst.push_back(QVector<int>());
        for (y=0;y<dest_y;y++)
        {
            fx=x*step_x;
            fy=y*step_y;
            dx=fx-(int)fx;
            dy=fy-(int)fy;
            getPixel(src,floor(fx),floor(fy),&R1);    // get N1 colours
            getPixel(src,ceil(fx),floor(fy),&R2); // get N2 colours
            getPixel(src,floor(fx),ceil(fy),&R3); // get N3 colours
            getPixel(src,ceil(fx),ceil(fy),&R4);  // get N4 colours
            // Interpolate to get T1 and T2 colours
            RT1=(dx*R2)+(1-dx)*R1;

            RT2=(dx*R4)+(1-dx)*R3;
            // Obtain final colour by interpolating between T1 and T2
            R=(unsigned char)((dy*RT2)+((1-dy)*RT1));

            // Store the final colour
            dst[x].push_back(R);
        }
    }
}


// Task 1:
void MainWindow::on_pushButton_2_clicked()
{
    // Reset stage
    ui->label_error->setText("");

    QRadioButton *r1 = ui->radioButton;
    QRadioButton *r2 = ui->radioButton_2;
    QRadioButton *r3 = ui->radioButton_3;
    QRadioButton *r4 = ui->radioButton_4;

    int destWidth = ui->lineEdit->text().toInt();
    int destHeight = ui->lineEdit_2->text().toInt();

    QVector<QVector<int>> newImg;
    int w = imgArray.size();
    int h = imgArray[0].size();

    if(r1->isChecked()) {
        // Linear x-axis
        ui->label_error->setText("Generating...");
        handleLinearInterpolationX(imgArray, newImg, w, h, destWidth, destHeight);
        ui->label_error->setText("");

    } else if(r2->isChecked()) {
        // Linear y-axis
        ui->label_error->setText("Generating...");
        handleLinearInterpolationY(imgArray, newImg, w, h, destWidth, destHeight);
        ui->label_error->setText("");

    } else if(r3->isChecked()) {
        // Bilinear
        ui->label_error->setText("Generating...");
        handleBilinearInterpolation(imgArray, newImg, w, h, destWidth, destHeight);
        ui->label_error->setText("");

    } else if(r4->isChecked()) {
        // Nearest Neighbors
        //Use nearest sampling

        double scale_w = w / (double)destWidth;
        double scale_h = h / (double)destHeight;

        int tSrcH = 0, tSrcW = 0;
        int index_src = 0, index_dest = 0;

        for (int y=0; y < destWidth; y++) {
            newImg.push_back(QVector<int>());

            for (int x=0; x < destHeight; x++) {
                newImg[y].push_back(imgArray[y* scale_w][x * scale_h]);
            }
        }

    }
    QPixmap result;
    grayscaleToQPixmap(newImg , result);
    ui->label_pic_2->setPixmap(result);
    if(destWidth < 100 && destHeight < 100) {
        ui->label_pic_3->setPixmap(result.scaled(ui->label_pic_3->width(),ui->label_pic_3->height()));
    }

    ui->label_pic->adjustSize();


    // Showing new size
    ui->label_new_w->setText(QString::number(newImg.size()));
    ui->label_new_h->setText(QString::number(newImg[0].size()));

    // Reset stage
    for(int i=0 ; i< newImg.size(); i++) {
        for(int j=0 ; j< newImg[i].size(); j++) {
            newImg[i].pop_back();
        }
        newImg.pop_front();
    }
    //    QString fileSave = "/Users/datnguyen/Downloads/ImgResult/Task1/result" + QString::number(rand() % 1000)+".png";
    //    result.save(fileSave, "PNG");
}

void handleReducingBit(QVector<QVector<int>> & src, QVector<QVector<int>> & newImg, int bits) {
    int num_colors = pow(2, bits);
    int divisor = 256 / num_colors;
    int max_quantized_value = 255 / divisor;

    for(int i=0 ; i< src.size(); i++) {
        newImg.push_back(QVector<int>());
        for(int j=0 ; j< src[i].size(); j++) {
            int new_value = ((imgArray[i][j] / divisor) * 255) / max_quantized_value;
            newImg[i].push_back(new_value);
        }
    }


}

// Task 2:
void MainWindow::on_pushButton_3_clicked()
{
    int newBit = ui->spinBox->text().toInt();
    QVector<QVector<int>> newImg;

    handleReducingBit(imgArray, newImg, newBit);



    // Showing result image
    int destWidth = ui->lineEdit->text().toInt();
    int destHeight = ui->lineEdit_2->text().toInt();
    QPixmap result;
    grayscaleToQPixmap(newImg , result);
    ui->label_pic_2->setPixmap(result);
    if(destWidth < 100 && destHeight < 100) {
        ui->label_pic_3->setPixmap(result.scaled(ui->label_pic_3->width(),ui->label_pic_3->height()));
    }

    ui->label_pic->adjustSize();
    // Showing new size
    ui->label_new_w->setText(QString::number(newImg.size()));
    ui->label_new_h->setText(QString::number(newImg[0].size()));

    // Reset stage
    for(int i=0 ; i< newImg.size(); i++) {
        for(int j=0 ; j< newImg[i].size(); j++) {
            newImg[i].pop_back();
        }
        newImg.pop_front();
    }
    //    imgArray = newImg;

    QString fileSave = "/Users/datnguyen/Downloads/ImgResult/Task2/result" + QString::number(rand() % 1000)+".png";
    result.save(fileSave, "PNG");
}


// Increase the zooming by 2
void MainWindow::on_pushButton_5_clicked()
{
    int w = ui->lineEdit->text().toInt();
    int h = ui->lineEdit_2->text().toInt();
    ui->lineEdit->setText(QString::number(w*2));
    ui->lineEdit_2->setText(QString::number(h*2));
}


void MainWindow::on_pushButton_4_clicked()
{
    int w = ui->lineEdit->text().toInt();
    int h = ui->lineEdit_2->text().toInt();
    if(w/2 != 0 && h/2 != 0) {
        ui->lineEdit->setText(QString::number(w/2));
        ui->lineEdit_2->setText(QString::number(h/2));
    } else {
        ui->label_error->setText("Cannot make smaller. Please try another dimensions");
    }
}




// ---------------------------------------------------------------- HOMEWORK 2
// Global histogram equalization
void MainWindow::on_pushButton_10_clicked()
{
    // Declaring 2 arrays for storing histogram values (frequencies) and
    // new gray level values (newly mapped pixel values as per algorithm)
    int hist[256] = { 0 };
    int new_gray_level[256] = { 0 };
    int input_file, output_file, col, row, total, curr, i;
    int rows = imgArray.length();
    int cols = imgArray[0].length();

    for (row = 0; row < rows; row++) {

        // logic for calculating histogram
        for (col = 0; col < cols; col++)
            hist[(int)imgArray[row][col]]++;
    }

    // calculating total number of pixels
    total = cols * rows;

    curr = 0;

    // calculating cumulative frequency and new gray levels
    for (i = 0; i < 256; i++) {
        // cumulative frequency
        curr += hist[i];

        // calculating new gray level after multiplying by
        // maximum gray count which is 255 and dividing by
        // total number of pixels
        new_gray_level[i] = round((((float)curr) * 255) / total);
    }

    QVector<QVector<int>> newImg;
    // performing histogram equalisation by mapping new gray levels
    for (row = 0; row < rows; row++) {
        // reading a row of image
        newImg.push_back(QVector<int> ());
        // mapping to new gray level values
        for (col = 0; col < cols; col++)
            newImg[row].push_back((unsigned char)new_gray_level[(int)imgArray[row][col]]);
    }

    handleRenderImage(newImg);
}

int part(const QVector<QVector<int>>& local) {
    QVector<int> possibility(256, 0);
    for (int v = 0; v < 3; v++) {
        for (int h = 0; h < 3; h++) {
            int P = local[v][h];
            possibility[P]++;
        }
    }
    for (int o = 0; o < 256; o++) {
        possibility[o] = static_cast<int>(std::round(static_cast<float>(possibility[o]) / 9 * 255));
    }
    int S = local[1][1];
    return possibility[S];
}


// Function to add duplicate value padding to a 2D array
QVector<QVector<int>> addDuplicatePadding(const QVector<QVector<int>>& inputImage, int paddingSize) {
    // Get the dimensions of the input image
    int height = inputImage.size();
    int width = inputImage[0].size();

    // Create a new 2D array with padding
    QVector<QVector<int>> paddedImage(height + 2 * paddingSize, QVector<int>(width + 2 * paddingSize, 0));

    // Copy the input image to the central region of the padded image
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            paddedImage[i + paddingSize][j + paddingSize] = inputImage[i][j];
        }
    }

    // Duplicate value padding for the top and bottom borders
    for (int i = 0; i < paddingSize; i++) {
        for (int j = 0; j < width + 2 * paddingSize; j++) {
            paddedImage[i][j] = paddedImage[paddingSize][j]; // Copy from the first row
            paddedImage[height + 2 * paddingSize - 1 - i][j] = paddedImage[height + paddingSize - 1][j]; // Copy from the last row
        }
    }

    // Duplicate value padding for the left and right borders
    for (int i = 0; i < height + 2 * paddingSize; i++) {
        for (int j = 0; j < paddingSize; j++) {
            paddedImage[i][j] = paddedImage[i][paddingSize]; // Copy from the first column
            paddedImage[i][width + 2 * paddingSize - 1 - j] = paddedImage[i][width + paddingSize - 1]; // Copy from the last column
        }
    }

    return paddedImage;
}


// Local histogram equalization
void MainWindow::on_pushButton_9_clicked()
{
    // Define the size of the filter
    int size = ui->lineEdit_3->text().toInt();

    ui->label_error->setText("Generating...");
    // Display the original histogram
    int heights = imgArray.length();
    int widths = imgArray[0].length();

    QVector<int> store(256, 0);
    QVector<int> saved(256, 0);
    for (int i = 0; i < heights; i++) {
        for (int j = 0; j < widths; j++) {
            int k = static_cast<int>(imgArray[i][j]);
            store[k]++;
        }
    }
    int padding = size / 2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(imgArray, padding);

    // Get the new image size
    int height_new = imgPad.length();
    int width_new = imgPad[0].length();
    std::cout << height_new << std::endl;
    std::cout << width_new << std::endl;

    for (int x = 0; x < 256; x++) {
        store[x] = 0;
        saved[x] = 0;
    }

    // Copy the original image, used to display the new result
    QVector<QVector<int>> imgArray_copy;

    // Looping for the entire image
    for (int a = padding; a <= height_new - padding -1; a++) {
        imgArray_copy.push_back(QVector<int>());
        for (int b = padding; b <= width_new - padding -1; b++) {
            // Looping for the entire window size
            for (int c = a - padding; c <= a + padding; c++) {
                for (int d = b - padding; d <= b + padding; d++) {
                    // Get the pixel value in the current pixel position
                    int k = static_cast<int>(imgPad[c][d]);
                    store[k]++;
                }
            }

            // Perform the cumulative distribution function
            QVector<int> sum_hist(256);
            sum_hist[0] = store[0];
            for (int x = 1; x < 256; x++) {
                sum_hist[x] = sum_hist[x - 1] + store[x];
            }

            // Get the new pixel value from the cumulative distribution
            for (int x = 0; x < 256; x++) {
                saved[x] = sum_hist[x] * 255 / sum_hist[255];
            }

            // Write the new pixel value into the copied image
            int k = static_cast<int>(imgPad[a][b]);
            imgArray_copy[a-padding].push_back(saved[k]);

            // Reset the stored value
            for (int x = 0; x < 256; x++) {
                store[x] = 0;
                saved[x] = 0;
            }
        }
    }

    handleRenderImage(imgArray_copy);
    ui->label_error->setText("");


}




QVector<QVector<int>> median_filter(QVector<QVector<int>> input, int filter_size) {
    QVector<QVector<int>> output;

    for (int i = 0; i < input.size(); i++) {
        output.push_back(QVector<int>());
        for (int j = 0; j < input[i].size(); j++) {
            QVector<int> window;
            for (int k = -filter_size / 2; k <= filter_size / 2; k++) {
                for (int l = -filter_size / 2; l <= filter_size / 2; l++) {
                    int ii = i + k;
                    int jj = j + l;
                    if (ii >= 0 && ii < input.size() && jj >= 0 && jj < input[ii].size()) {
                        window.push_back(input[ii][jj]);
                    }
                }
            }

            sort(window.begin(), window.end());
            output[i].push_back(window[window.size() / 2]);
            window.clear();
        }
    }

    return output;
}


QVector<QVector<int>> filter2dArrayImage(const QVector<QVector<int>>& image, QVector<QVector<double>>  filter) {


    // Get the dimensions of the filter.
    int filterRows = filter.length();
    int filterColumns = filter[0].length();

    int padding = filterRows/2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(image, padding);


    QVector<QVector<int>> output;
    // Iterate over the pixels of the input image.
    for (int a = 0; a < image.length(); a++) {
        output.push_back(QVector<int>());
        for (int b = 0; b < image[0].length(); b++) {
            // Looping for the entire window size
            int filteredPixelValue = 0;

            for (int c = 0; c < filterRows; c++) {
                for (int d = 0; d < filterColumns; d++) {
                    // Get the pixel value in the current pixel position
                    filteredPixelValue += (int)(filter[c][d] * imgPad[a+c][b+d]);
                }
            }
            output[a].push_back(filteredPixelValue);
        }
    }
    return output;
}

// Gaussian Smoothing Filter
QVector<QVector<double>> smooth_filter(int size) {
    QVector<QVector<double>> kernel;
    double sigma = 2.0;
    double p, q = 2.0 * sigma * sigma;
    double sum = 0.0;
    int d = size/2;
    double PI  = 3.1415;
    for (int x = -d; x <= d; x++) {
        kernel.push_back(QVector<double>());
        for (int y = -d; y <= d; y++) {
            p = sqrt(x * x + y * y);
            kernel[x + d].push_back(exp((-(p * p) / q)) / (PI * q));
            sum += kernel[x + d][y + d];
        }
    }
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            kernel[i][j] /= sum;
    return kernel;
}



// All filter handling
void MainWindow::on_pushButton_8_clicked()
{
    ui->label_error->setText("");
    int size = ui->lineEdit_3->text().toInt();

    QRadioButton *r1 = ui->radioButton_5;
    QRadioButton *r2 = ui->radioButton_6;
    QRadioButton *r3 = ui->radioButton_7;
    QRadioButton *r4 = ui->radioButton_8;

    QVector<QVector<int>> newImg;

    ui->label_error->setText("Generating...");

    if(r1->isChecked()) {
        QVector<QVector<double>> smoothFil = smooth_filter(size);
        newImg = filter2dArrayImage(imgArray, smoothFil);

    } else if(r2->isChecked()) {
        newImg = median_filter(imgArray, size);

    } else if(r3->isChecked()) {
        // Laplacian Sharpening Filter (8-neighbors)
        QVector<QVector<int>> laplacianFilter = {{0,1,0}, {1,-4,1}, {0,1,0}};
        for (int y = 0; y < imgArray.length(); y++) {
            newImg.push_back(QVector<int>());
            for (int x = 0; x < imgArray[0].length(); x++) {
                newImg[y].push_back(0);
            }
        }
        int height = imgArray.length();
        int width = imgArray[0].length();
        int filterHeight = 3;
        int filterWidth = 3;
        int newImageHeight = height - filterHeight + 1;
        int newImageWidth = width - filterWidth + 1;
        int i, j, h, w;

        //convolution
        for (i = 0; i < newImageHeight; i++) {
            for (j = 0; j < newImageWidth; j++) {
                int temp = newImg[i][j];
                for (h = i; h < i + filterHeight; h++) {
                    for (w = j; w < j + filterWidth; w++) {
                        temp += laplacianFilter[h - i][w - j] * imgArray[h][w];
                    }
                }
                newImg[i][j]= temp;
            }
        }

        //img - laplace
        for (int y = 0; y < newImg.length(); y++) {
            for (int x = 0; x < newImg[y].length(); x++) {
                newImg[y][x] = (imgArray[y][x] - newImg[y][x]);
            }
        }

    } else if(r4->isChecked()) {
        for (int y = 0; y < imgArray.length(); y++) {
            newImg.push_back(QVector<int>());
            for (int x = 0; x < imgArray[0].length(); x++) {
                newImg[y].push_back(0);
            }
        }
        QVector<QVector<double>> smoothFil = smooth_filter(size);
        QVector<QVector<int>> blurImg = filter2dArrayImage(imgArray, smoothFil);
        int k = ui->lineEdit_5->text().toInt();
        for (int y = 0; y < newImg.length(); y++) {
            for (int x = 0; x < newImg[y].length(); x++) {
                newImg[y][x] = imgArray[y][x] + k*(imgArray[y][x] - blurImg[y][x]);
            }
        }


    }

    handleRenderImage(newImg);
    ui->label_error->setText("");

}
// Bit plane
void slicingBitPlane(int d, QVector<QVector<int>> & newImg) {
    int bitsPerPixel = 8;

    for (int i = 0; i < imgArray.length(); i++)
        for (int j = 0; j < imgArray[i].length(); j++) {

            int col = 0;
            int fcol = imgArray[i][j];
            if (((fcol >> d) & 1) > 0) col = 1;
            newImg[i][j] = 255*col;
        }
}

// Render bit plane

void MainWindow::on_pushButton_11_clicked()
{
    int newBit = ui->spinBox_2->text().toInt();
    QVector<QVector<int>> newImg;
    for (int y = 0; y < imgArray.length(); y++) {
        newImg.push_back(QVector<int>());
        for (int x = 0; x < imgArray[0].length(); x++) {
            newImg[y].push_back(0);
        }
    }
    slicingBitPlane(newBit, newImg);
    handleRenderImage(newImg);


//    QPixmap result;
//    grayscaleToQPixmap(newImg , result);
//    QString fileSave = "/Users/datnguyen/Downloads/ImgResult/result" + QString::number(newBit)+".png";
//    result.save(fileSave, "PNG");
}

// Recombine image from lowest bit plane
void MainWindow::on_pushButton_12_clicked()
{
    int d = ui->spinBox_3->text().toInt();
    QVector<QVector<int>> plane;
    for (int y = 0; y < imgArray.length(); y++) {
        plane.push_back(QVector<int>());
        for (int x = 0; x < imgArray[0].length(); x++) {
            plane[y].push_back(0);
        }
    }
    QVector<QVector<int>> newImg;
    for (int y = 0; y < imgArray.length(); y++) {
        newImg.push_back(QVector<int>());
        for (int x = 0; x < imgArray[0].length(); x++) {
            newImg[y].push_back(0);
        }
    }
    for(int k = 0; k <d; k++) {
        for (int i = 0; i < imgArray.length(); i++){
            for (int j = 0; j < imgArray[i].length(); j++) {

                int col = 0;
                int fcol = imgArray[i][j];
                if (((fcol >> k) & 1) > 0) col = 1;
                plane[i][j] = col;
            }
        }
        for (int i = 0; i < imgArray.length(); i++)
            for (int j = 0; j < imgArray[i].length(); j++) {
                newImg[i][j] += (plane[i][j] * pow(2,k));
                plane[i][j]=0;
            }
    }

    handleRenderImage(newImg);
}


void MainWindow::on_pushButton_13_clicked()
{
    int d = ui->spinBox_3->text().toInt();
    QVector<QVector<int>> plane;
    for (int y = 0; y < imgArray.length(); y++) {
        plane.push_back(QVector<int>());
        for (int x = 0; x < imgArray[0].length(); x++) {
            plane[y].push_back(0);
        }
    }
    QVector<QVector<int>> newImg;
    for (int y = 0; y < imgArray.length(); y++) {
        newImg.push_back(QVector<int>());
        for (int x = 0; x < imgArray[0].length(); x++) {
            newImg[y].push_back(0);
        }
    }
    for(int k = 7; k >=8-d; k--) {
        for (int i = 0; i < imgArray.length(); i++){
            for (int j = 0; j < imgArray[i].length(); j++) {

                int col = 0;
                int fcol = imgArray[i][j];
                if (((fcol >> k) & 1) > 0) col = 1;
                plane[i][j] = col;
            }
        }
        for (int i = 0; i < imgArray.length(); i++)
            for (int j = 0; j < imgArray[i].length(); j++) {
                newImg[i][j] += (plane[i][j] * pow(2,k));
                plane[i][j]=0;
            }
    }

    handleRenderImage(newImg);
}




// Change filter size

void MainWindow::on_pushButton_7_clicked()
{
    int w = ui->lineEdit_3->text().toInt();
    int h = ui->lineEdit_4->text().toInt();
    ui->lineEdit_3->setText(QString::number(w+2));
    ui->lineEdit_4->setText(QString::number(h+2));
}


void MainWindow::on_pushButton_6_clicked()
{
    int w = ui->lineEdit_3->text().toInt();
    int h = ui->lineEdit_4->text().toInt();
    if(w-2 != 0 && h-2 != 0) {
        ui->lineEdit_3->setText(QString::number(w-2));
        ui->lineEdit_4->setText(QString::number(h-2));
    } else {
        ui->label_error->setText("Cannot make smaller. Please try another dimensions");
    }
}


void MainWindow::on_lineEdit_3_textChanged(const QString &arg1)
{
    ui->lineEdit_4->setText(arg1);
}



void MainWindow::on_lineEdit_4_textEdited(const QString &arg1)
{
    ui->lineEdit_3->setText(arg1);
}







// -------------------------------------------------------------------- Assignment 3 ---------------------------------------------------------

QVector<QVector<int>> contraharmonicFilter2dArrayImage(const QVector<QVector<int>>& image, int filterRows, int filterColumns, int Q) {

    int padding = filterRows/2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(image, padding);


    QVector<QVector<int>> output;
    // Iterate over the pixels of the input image.
    for (int a = 0; a < image.length(); a++) {
        output.push_back(QVector<int>());
        for (int b = 0; b < image[0].length(); b++) {
            // Looping for the entire window size
            double filteredPixelValue1 = 0;
            double filteredPixelValue2 = 0;

            for (int c = 0; c < filterRows; c++) {
                for (int d = 0; d < filterColumns; d++) {
                    // Get the pixel value in the current pixel position
                    filteredPixelValue1 += (double)pow(imgPad[a+c][b+d], Q);
                    filteredPixelValue2 += (double)pow(imgPad[a+c][b+d], Q+1);
                }
            }
            output[a].push_back((int)filteredPixelValue2/filteredPixelValue1);
        }
    }
    return output;
}


QVector<QVector<int>> geometricFilter2dArrayImage(const QVector<QVector<int>>& image, int filterRows, int filterColumns) {

    int padding = filterRows/2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(image, padding);


    QVector<QVector<int>> output;
    // Iterate over the pixels of the input image.
    for (int a = 0; a < image.length(); a++) {
        output.push_back(QVector<int>());
        for (int b = 0; b < image[0].length(); b++) {
            // Looping for the entire window size
            double filteredPixelValue = 0;

            for (int c = 0; c < filterRows; c++) {
                for (int d = 0; d < filterColumns; d++) {
                    // Get the pixel value in the current pixel position
                    filteredPixelValue += log(imgPad[a+c][b+d]);
                }
            }
            output[a].push_back(exp(filteredPixelValue/ (filterRows*filterColumns)));
        }
    }
    return output;
}

QVector<QVector<int>> maxFilter2dArrayImage(const QVector<QVector<int>>& image, int filterRows, int filterColumns) {

    int padding = filterRows/2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(image, padding);


    QVector<QVector<int>> output;
    // Iterate over the pixels of the input image.
    for (int a = 0; a < image.length(); a++) {
        output.push_back(QVector<int>());
        for (int b = 0; b < image[0].length(); b++) {
            // Looping for the entire window size
            int maxx = 0;

            for (int c = 0; c < filterRows; c++) {
                for (int d = 0; d < filterColumns; d++) {
                    // Get the pixel value in the current pixel position
                    if(imgPad[a+c][b+d] >= maxx) {
                        maxx = imgPad[a+c][b+d];
                    }
                }
            }
            output[a].push_back(maxx);
        }
    }
    return output;
}

QVector<QVector<int>> minFilter2dArrayImage(const QVector<QVector<int>>& image, int filterRows, int filterColumns) {

    int padding = filterRows/2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(image, padding);


    QVector<QVector<int>> output;
    // Iterate over the pixels of the input image.
    for (int a = 0; a < image.length(); a++) {
        output.push_back(QVector<int>());
        for (int b = 0; b < image[0].length(); b++) {
            // Looping for the entire window size
            int minn = 256;

            for (int c = 0; c < filterRows; c++) {
                for (int d = 0; d < filterColumns; d++) {
                    // Get the pixel value in the current pixel position
                    if(imgPad[a+c][b+d] <= minn) {
                        minn = imgPad[a+c][b+d];
                    }
                }
            }
            output[a].push_back(minn);
        }
    }
    return output;
}

QVector<QVector<int>> midPointFilter2dArrayImage(const QVector<QVector<int>>& image, int filterRows, int filterColumns) {

    int padding = filterRows/2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(image, padding);


    QVector<QVector<int>> output;
    // Iterate over the pixels of the input image.
    for (int a = 0; a < image.length(); a++) {
        output.push_back(QVector<int>());
        for (int b = 0; b < image[0].length(); b++) {
            // Looping for the entire window size
            int maxx = 0;
            int minn = 256;

            for (int c = 0; c < filterRows; c++) {
                for (int d = 0; d < filterColumns; d++) {
                    // Get the pixel value in the current pixel position
                    if(imgPad[a+c][b+d] >= maxx) {
                        maxx = imgPad[a+c][b+d];
                    }
                    if(imgPad[a+c][b+d] <= minn) {
                        minn = imgPad[a+c][b+d];
                    }
                }
            }
            output[a].push_back((int) (0.5 * (maxx + minn)));
        }
    }
    return output;
}

QVector<QVector<int>> alphaTrimmedFilter2dArrayImage(const QVector<QVector<int>>& image, int filterRows, int filterColumns, int d) {

    int padding = filterRows/2;
    QVector<QVector<int>> imgPad = addDuplicatePadding(image, padding);


    QVector<QVector<int>> output;
    // Iterate over the pixels of the input image.
    for (int a = 0; a < image.length(); a++) {
        output.push_back(QVector<int>());
        for (int b = 0; b < image[0].length(); b++) {
            // Looping for the entire window size
            QVector<int> tmp;
            for (int c = 0; c < filterRows; c++) {
                for (int d = 0; d < filterColumns; d++) {
                    tmp.push_back( imgPad[a+c][b+d]);
                }
            }
            sort(tmp.begin(), tmp.end());
            int total = 0;
            for (int i = 0; i < tmp.size(); i++) {
                if(i > d/2 && i < tmp.size() - 1 - d/2) {
                    total += tmp[i];
                }
            }
            output[a].push_back((int) total/(filterRows*filterColumns - d));
            tmp.clear();
        }
    }
    return output;
}
void MainWindow::on_pushButton_16_clicked()
{
    ui->label_error->setText("");

    QRadioButton *r1 = ui->radioButton_9;
    QRadioButton *r2 = ui->radioButton_10;
    QRadioButton *r3 = ui->radioButton_11;
    QRadioButton *r4 = ui->radioButton_12;
    QRadioButton *r5 = ui->radioButton_13;
    QRadioButton *r6 = ui->radioButton_14;
    QRadioButton *r7 = ui->radioButton_15;
    QRadioButton *r8 = ui->radioButton_16;


    int filterW = ui->lineEdit_6->text().toInt();
    int filterH = ui->lineEdit_7->text().toInt();

    QVector<QVector<int>> newImg;
    ui->label_error->setText("Generating...");

    if(r1->isChecked()) {
        newImg = contraharmonicFilter2dArrayImage(imgArray, filterW, filterH, 0);
    } else if(r2->isChecked()) {
        newImg = geometricFilter2dArrayImage(imgArray, filterW, filterH);
    } else if(r3->isChecked()) {
        newImg = contraharmonicFilter2dArrayImage(imgArray, filterW, filterH, -1);

    } else if(r4->isChecked()) {
        int Q = ui->lineEdit_Q->text().toInt();
        newImg = contraharmonicFilter2dArrayImage(imgArray, filterW, filterH, Q);

    }else if(r5->isChecked()) {
        newImg = maxFilter2dArrayImage(imgArray, filterW, filterH);
    }else if(r6->isChecked()) {
        newImg = minFilter2dArrayImage(imgArray, filterW, filterH);
    }else if(r7->isChecked()) {
        newImg = midPointFilter2dArrayImage(imgArray, filterW, filterH);
    }else if(r8->isChecked()) {
        int d = ui->lineEdit_d->text().toInt();
        newImg = alphaTrimmedFilter2dArrayImage(imgArray, filterW, filterH, d);
    }

    ui->label_error->setText("");

    QPixmap result;
    grayscaleToQPixmap(newImg , result);
    ui->label_pic_2->setPixmap(result);
    if(newImg.size() < 100 && newImg[0].size() < 100) {
        ui->label_pic_3->setPixmap(result.scaled(ui->label_pic_3->width(),ui->label_pic_3->height()));
    }

    ui->label_pic->adjustSize();


    // Showing new size
    ui->label_new_w->setText(QString::number(newImg.size()));
    ui->label_new_h->setText(QString::number(newImg[0].size()));

    // Reset stage
    for(int i=0 ; i< newImg.size(); i++) {
        for(int j=0 ; j< newImg[i].size(); j++) {
            newImg[i].pop_back();
        }
        newImg.pop_front();
    }
    //    QString fileSave = "/Users/datnguyen/Downloads/ImgResult/Task1/result" + QString::number(rand() % 1000)+".png";
    //    result.save(fileSave, "PNG");
}






// Change filter size
void MainWindow::on_pushButton_15_clicked()
{
    int w = ui->lineEdit_6->text().toInt();
    int h = ui->lineEdit_7->text().toInt();
    if(w-2 != 0 && h-2 != 0) {
        ui->lineEdit_6->setText(QString::number(w-2));
        ui->lineEdit_7->setText(QString::number(h-2));
    } else {
        ui->label_error->setText("Cannot make smaller. Please try another dimensions");
    }
}
void MainWindow::on_pushButton_14_clicked()
{
    int w = ui->lineEdit_6->text().toInt();
    int h = ui->lineEdit_7->text().toInt();
    ui->lineEdit_6->setText(QString::number(w+2));
    ui->lineEdit_7->setText(QString::number(h+2));
}


void MainWindow::on_lineEdit_6_textChanged(const QString &arg1)
{
    ui->lineEdit_7->setText(arg1);
}


void MainWindow::on_lineEdit_7_textChanged(const QString &arg1)
{
    ui->lineEdit_6->setText(arg1);
}


// -------------------------------------------------------------------- Assignment 4 ---------------------------------------------------------

QVector<QVector<QPair<int, int>>> RLEncode(const QVector<QVector<int>>& image) {
    int height = image.length();
    int width = image[0].length();

    // Create a vector to store the encoded image
    QVector<QVector<QPair<int, int>>> encodedImage;

    // Iterate over the image
    for (int i = 0; i < height; i++) {
        // Get the current pixel value
        int currentPixel = image[i][0];

        // Initialize the run length
        int runLength = 1;

        // Iterate over the remaining pixels in the row
        for (int j = 1; j < width; j++) {
            // If the current pixel value is the same as the previous pixel value, increment the run length
            if (currentPixel == image[i][j]) {
                runLength++;
            } else {
                // Otherwise, encode the run length and reset it
                encodedImage[i].push_back(qMakePair(currentPixel, runLength));
                currentPixel = image[i][j];
                runLength = 1;
            }
        }
        // Encode the last run length
        encodedImage[i].push_back(qMakePair(currentPixel, runLength));
    }

    return encodedImage;
}

QVector<QVector<int>> RLDecode(const QVector<QVector<QPair<int, int>>> encodedImage) {
    // Get the dimensions of the encoded image
    int height = encodedImage.length();

    // Create a vector to store the decoded image
    QVector<QVector<int>> decodedImage;

    // Iterate over the encoded image
    for (int i = 0; i < height; i++) {
        decodedImage.push_back(QVector<int>());

        for (int j = 0; j < encodedImage[i].length(); j++) {
            int runLength = encodedImage[i][j].second;

            // Initialize the current pixel value
            int currentPixel = encodedImage[i][j].first;

            // Iterate over the remaining pixels in the row
            for (int j = 0; j < runLength; j++) {
                // Decode the run length and reset it
                decodedImage[i].push_back(currentPixel);
                runLength--;
            }

        }
    }

    return decodedImage;
}

QVector<QVector<int>> runLengthCodingOnGrayscale(const QVector<QVector<int>>& image) {
    QVector<QVector<QPair<int, int>>> encodedImage = RLEncode(image);
    for (int i = 0; i < encodedImage.length(); i++) {

        for (int j = 0; j < encodedImage[i].length(); j++) {
            int runLength = encodedImage[i][j].second;
            int currentPixel = encodedImage[i][j].first;
            cout << runLength << " " << currentPixel << " ";
        }
        cout << endl;
    }
    return QVector<QVector<int>>();
//    return RLDecode(encodedImage);
}

//QVector<QVector<int>> runLengthCodingOnBitPlanes(const QVector<QVector<int>>& image) {

//}

//QVector<QVector<int>> variableLengthHuffmanCoding(const QVector<QVector<int>>& image) {

//}
//QVector<QVector<int>> LZW(const QVector<QVector<int>>& image) {

//}

void MainWindow::on_pushButton_17_clicked()
{
    ui->label_error->setText("");

    QRadioButton *r1 = ui->radioButton_17;
    QRadioButton *r2 = ui->radioButton_18;
    QRadioButton *r3 = ui->radioButton_19;
    QRadioButton *r4 = ui->radioButton_20;

    QVector<QVector<int>> newImg;
    ui->label_error->setText("Generating...");

    if(r1->isChecked()) {
        newImg = runLengthCodingOnGrayscale(imgArray);
    }
//    else if(r2->isChecked()) {
//        newImg = runLengthCodingOnBitPlanes(imgArray);

//    } else if(r3->isChecked()) {
//        newImg = variableLengthHuffmanCoding(imgArray);

//    } else if(r4->isChecked()) {
//        newImg = LZW(imgArray);
//    }

//    ui->label_error->setText("");

//    QPixmap result;
//    grayscaleToQPixmap(newImg , result);
//    ui->label_pic_2->setPixmap(result);
//    if(newImg.size() < 100 && newImg[0].size() < 100) {
//        ui->label_pic_3->setPixmap(result.scaled(ui->label_pic_3->width(),ui->label_pic_3->height()));
//    }

//    ui->label_pic->adjustSize();


//    // Showing new size
//    ui->label_new_w->setText(QString::number(newImg.size()));
//    ui->label_new_h->setText(QString::number(newImg[0].size()));

//    // Reset stage
//    for(int i=0 ; i< newImg.size(); i++) {
//        for(int j=0 ; j< newImg[i].size(); j++) {
//            newImg[i].pop_back();
//        }
//        newImg.pop_front();
//    }
}

