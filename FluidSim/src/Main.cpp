#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <vector>
#include <random>
#include "stb/stb_image.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);

int SCREEN_WIDTH = 256;
int SCREEN_HEIGHT = 256;

float rectX = -1.0f;
float rectY = -1.0f;
float rectWidth = 1.0f;
float rectHeight = 1.0f;
bool isDragging = false;
float previousX = 0.0;
float previousY = 0.0;

int N = 64;
int iter = 4;

int temps = 0;
int t = 0;

static const char* vertexShaderSource = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 2) in vec2 aTexCoord;

out vec2 TexCoord;
void main()
{
gl_Position = vec4(aPos, 1.0);
TexCoord = aTexCoord;
}
    )";

static const char* fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;

in vec2 TexCoord;
uniform sampler2D densityTexture;
void main()
{
vec4 texColor = texture(densityTexture, TexCoord);
FragColor = texColor;
//FragColor = vec4(texColor.rgb, 1.0);
}
    )";

int IX(int x, int y) {
    int xlimit = std::min(std::max(x, 0), N-1);
    int ylimit = std::min(std::max(y, 0), N-1);
    int index = (xlimit + N * ylimit);
    return index;
};

float clamp(float input, float lower, float upper) {
    return std::min(std::max(input, lower), upper - 1);
}

struct fluidGridType {

    int N;
    float visc;
    float dt;
    float diff;

    std::vector<float> Vx;
    std::vector<float> Vy;

    std::vector<float> Vx0;
    std::vector<float> Vy0;

    std::vector<float> density;
    std::vector<float> s;
};

//fluidGridType fluidGridCreate(fluidGridType &newGrid) {
fluidGridType fluidGridCreate() {
    fluidGridType newGrid;
    newGrid.N = N;
    newGrid.dt = 0.2;
    newGrid.diff = 0.0;
    newGrid.visc = 0.0000001;
    int gridSize = newGrid.N * newGrid.N;
    newGrid.Vx.assign(gridSize, 0.0);
    newGrid.Vy.assign(gridSize, 0.0);
    newGrid.Vx0.assign(gridSize, 0.0);
    newGrid.Vy0.assign(gridSize, 0.0);
    newGrid.density.assign(gridSize, 0.0);
    newGrid.s.assign(gridSize, 0.0);
    return newGrid;
};

void addDensity(fluidGridType* grid, int x, int y, float amount) {
    grid->density[IX(x, y)] += amount;
    //grid->density[IX(x, y)] += grid->dt*grid->s[IX(x, y)];
};

void addVelocity(fluidGridType* grid, int x, int y, float amountX, float amountY) {
    int index = IX(x, y);
    grid->Vx[index] += amountX;
    grid->Vy[index] += amountY;

    //grid->Vx[index] += grid->dt * grid->Vx0[IX(x, y)];
    //grid->Vy[index] += grid->dt * grid->Vy0[IX(x, y)];
};

static void set_bnd(int b, std::vector<float>& x) {

    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }

    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)]
        + x[IX(0, 1)]);

    x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)]
        + x[IX(0, N - 2)]);

    x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)]
        + x[IX(N - 1, 1)]);

    x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)]
        + x[IX(N - 1, N - 2)]);
};

static void lin_solve(int b, std::vector<float> &x, std::vector<float> &x0, float a, float c) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip;
            }
        }
        set_bnd(b, x);
    }
};

static void diffuse(int b, std::vector<float> &x, std::vector<float>& x0, float diff, float dt) {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 4 * a);
};

static void project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div) {
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f * (
                velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] + velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]
                ) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 4);

    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= -0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
        }
    }

    set_bnd(1, velocX);
    set_bnd(2, velocY);
};

static void advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY, float dt) {
    float i0, i1, j0, j1;
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;

    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];

            x = i - dt * velocX[IX(i, j)];
            y = j - dt * velocY[IX(i, j)];

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
            i0 = int(x);
            i1 = i0 + 1.0f;
            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
            j0 = int(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;


            d[IX(i, j)] =
                s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }

    set_bnd(b, d);
};

static void step(fluidGridType* grid) {
    diffuse(1, grid->Vx0, grid->Vx, grid->visc, grid->dt);
    diffuse(2, grid->Vy0, grid->Vy, grid->visc, grid->dt);

    project(grid->Vx0, grid->Vy0, grid->Vx, grid->Vy);

    advect(1, grid->Vx, grid->Vx0, grid->Vx0, grid->Vy0, grid->dt);
    advect(2, grid->Vy, grid->Vy0, grid->Vx0, grid->Vy0, grid->dt);

    project(grid->Vx, grid->Vy, grid->Vx0, grid->Vy0);

    diffuse(0, grid->s, grid->density, grid->diff, grid->dt);
    advect(0, grid->density, grid->s, grid->Vx, grid->Vy, grid->dt);
};

static void test(fluidGridType* grid) {
    int min = 50;
    int max = 100;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(min, max);
    int amtRand = distrib(gen);
    rectX = distrib(gen);
    rectY = distrib(gen);
    float velX = float(distrib(gen));
    float velY = float(distrib(gen));
    //addDensity(grid, rectX, rectY, 150.0);
    //addVelocity(grid, velX, velY, 0.1, 0.1);

    float cx = float(0.5 * N);
    float cy = float(0.5 * N);

    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            addDensity(grid, cx + i, cy + j, amtRand);
        }
    }
    for (int i = 0; i < 2; i++) {
        static std::default_random_engine generator;
        static std::uniform_real_distribution<float> distribution(0.0, 1.0);
        float noise = distribution(generator);
        float angle = noise * 6.28318530718 * 2;
        float vx = cos(angle) * 0.2;
        float vy = sin(angle) * 0.2;
        addVelocity(grid, cx, cy, 0.1, 0.1);
    }
}


static void draw(fluidGridType* grid, unsigned int texture) {

    test(grid);
    step(grid);

    std::vector<unsigned char> textureData(N * N * 4);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float d = grid->density[IX(i, j)];
            textureData[4 * IX(i, j)] = d;  // R
            textureData[4 * IX(i, j) + 1] = d;  // G
            textureData[4 * IX(i, j) + 2] = d;  // B
            textureData[4 * IX(i, j) + 3] = 255.0;
        }
    }

    glBindTexture(GL_TEXTURE_2D, texture);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, N, N, GL_RGBA, GL_UNSIGNED_BYTE, textureData.data());
    //glBindTexture(GL_TEXTURE_2D, texture);
}


int main(void)

{

    fluidGridType fluidGrid = fluidGridCreate();

    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    /*
    glfwCreateWindow takes in width, height, name
        returns a GLFWwindow object
    */
    window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Hello World Triangle", NULL, NULL);
    if (!window)
    {
        std::cout << "Failed to create GLFW Window" << std::endl;
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    GLenum err = glewInit();
    if (err != GLEW_OK) {
        std::cout << "Error initializing GLEW" << std::endl;
        return -1;
    }

    glfwSetWindowUserPointer(window, &fluidGrid);

    /* We register callback functions after creating the window, but before the render loop */
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);

    int success;
    char infoLog[512];

    unsigned int vertexShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" <<
            infoLog << std::endl;
    }

    unsigned int fragmentShader;
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    //int success;
    //char infoLog[512];
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" <<
            infoLog << std::endl;
    }

    unsigned int shaderProgram;
    shaderProgram = glCreateProgram();
   
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::COMPILATION_FAILED\n" <<
            infoLog << std::endl;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    float vertices[] = {
        // positions // colors // texture coords
        0.5f, 0.5f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, // top right
        0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, // bottom right
        -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, // bottom left
        -0.5f, 0.5f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f // top left
    };

    unsigned int indices[] = {
    0, 1, 3, // first triangle
    1, 2, 3  // second triangle
    };

    unsigned int VBO, VAO, EBO;
    /* A vertex array object (VAO) can be bound like a vertex buffer object. It has the advantage of storting vertex attribute pointers (such as vertex shaders) so we only have to make the calls once */
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);


    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*) (6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    stbi_set_flip_vertically_on_load(true);
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    const void* dataPtr = static_cast<const void*>(fluidGrid.density.data());

    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, N, N, 0, GL_RGBA, GL_UNSIGNED_BYTE, dataPtr);

    // set the texture wrapping/filtering options (on currently bound texture)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // load and generate the texture
    int width, height, nrChannels;
    unsigned char* data = stbi_load("C:/Users/JTSte/Downloads/white_img.jpg", &width, &height,
        &nrChannels, 0);
    if (data)
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, N, N, 0, GL_RGBA,
            GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    }

    stbi_image_free(data);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {

        processInput(window);

        draw(&fluidGrid, texture);

        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);

        glActiveTexture(GL_TEXTURE0);

        glBindTexture(GL_TEXTURE_2D, texture);


        glUseProgram(shaderProgram);

        glBindVertexArray(VAO);

        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events (such as inputs) */
        glfwPollEvents();
    }

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);

    /* Clean all of GLFW's resources */
    glfwTerminate();
    return 0;
}
//
void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}
//
///* Resizes the viewpoint when the window is resized */
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
//
//    /* Tells OpenGL the size of the rendering window so it can display data & coords w.r.t the window */
    glViewport(0, 0, width, height);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if (action == GLFW_PRESS)
        {
            isDragging = true;
        }
        else if (action == GLFW_RELEASE)
        {
            isDragging = false;
        }
    }
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (isDragging)
    {
        // Convert mouse position to OpenGL coordinates
        rectX = (float)xpos;
        rectY = SCREEN_HEIGHT - (float)ypos; // Invert Y axis to match OpenGL coordinates

        rectX = clamp(rectX - 96, 0, N-1);
        rectY = clamp(rectY - 96, 0, N-1);

        float amtX = rectX - previousX;
        float amtY = rectY - previousY;

        fluidGridType* fluidGrid = static_cast<fluidGridType*>(glfwGetWindowUserPointer(window));
        addDensity(fluidGrid, rectX, rectY, 150.0);
        addVelocity(fluidGrid, rectX, rectY, amtX, amtY);
        previousX = rectX;
        previousY = rectY;
    }
}