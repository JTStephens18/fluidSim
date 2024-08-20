#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <iostream>
#include <vector>
#include "stb/stb_image.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);

int SCREEN_WIDTH = 64;
int SCREEN_HEIGHT = 64;

float rectX = -1.0f;
float rectY = -1.0f;
float rectWidth = 1.0f;
float rectHeight = 1.0f;
bool isDragging = false;

float verticess[] = {
    0.0f, 0.0f, // Bottom-left corner
    1.0f, 0.0f, // Bottom-right corner
    1.0f, 1.0f, // Top-right corner
    0.0f, 1.0f  // Top-left corner
};

// Indices for the rectangle (two triangles)
unsigned int indicess[] = {
    0, 1, 2,   // First triangle
    2, 3, 0    // Second triangle
};

static const char* vertexShaderSource1 = "#version 460 core\n"
"layout (location = 0) in vec3 aPos; \n"
"void main()\n"
"{\n"
"gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"}\0";

static const char* vertexShaderSources = R"(
        #version 330 core
        layout (location = 0) in vec2 aPos;

        uniform vec2 uScreenSize;
        uniform vec2 uPosition;
        uniform vec2 uSize;

        void main() {
            vec2 pos = aPos * uSize + uPosition;
            vec2 ndcPos = (pos / uScreenSize) * 2.0 - 1.0;
            gl_Position = vec4(ndcPos, 0.0, 1.0);
        }
    )";


static const char* vertexShaderSource_1 = R"(
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec3 aColor;

out vec3 ourColor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    gl_Position = projection * view * model * vec4(aPos, 1.0, 1.0);
    ourColor = aColor;
}
    )";

static const char* vertexShaderSource_2 = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoord;
out vec3 ourColor;
out vec2 TexCoord;
void main()
{
gl_Position = vec4(aPos, 1.0);
ourColor = aColor;
TexCoord = aTexCoord;
}
    )";

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

static const char* fragmentShaderSource1 = "#version 460 core\n"
"out vec4 FragColor;\n"
"void main()\n"
"{\n"
"FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);\n"
"}\0";

static const char* fragmentShaderSource2 = "#version 460 core\n"
"out vec4 FragColor;\n"
"void main()\n"
"{\n"
"FragColor = vec4(1.0f, 1.0f, 0.0f, 1.0f);\n"
"}\0";

static const char* fragmentShaderSource_1 = R"(
        #version 330 core
        
        uniform float uOpacity;
        uniform sampler2D gridTexture;

        out vec4 FragColor;
       
        void main() {
            FragColor = vec4(1.0, 1.0, 1.0, uOpacity); // White color
        }
    )";

static const char* fragmentShaderSource_2 = R"(
#version 330 core
out vec4 FragColor;
in vec3 ourColor;
in vec2 TexCoord;
uniform sampler2D ourTexture;
void main()
{
FragColor = texture(ourTexture, TexCoord);
}
    )";


static const char* fragmentShaderSource_3 = R"(
#version 330 core
out vec4 FragColor;
in vec3 ourColor;
in vec2 TexCoord;
uniform sampler2D densityTexture;
void main()
{
//float density = texture(densityTexture, TexCoord);
//FragColor = vec4(1.0, 1.0, 1.0, density);
FragColor = texture(densityTexture, TexCoord);
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

int N = 64;
int iter = 4;

int IX(int x, int y) {
    int xLimit = std::min(std::max(x, 0), N-1);
    int yLimit = std::min(std::max(y, 0), N-1);
    int index = xLimit + yLimit * N;
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

fluidGridType fluidGridCreate() {
    fluidGridType newGrid;
    newGrid.N = N;
    newGrid.dt = 0.1;
    newGrid.diff = 0;
    newGrid.visc = 0;
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
};

void addVelocity(fluidGridType* grid, int x, int y, float amountX, float amountY) {
    int index = IX(x, y);
    grid->Vx[index] += amountX;
    grid->Vy[index] += amountY;
};

static void set_bnd(int b, std::vector<float> x) {
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }

    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    x[IX(0, 0)] = 0.33f * (x[IX(1, 0)]
        + x[IX(0, 1)]
        + x[IX(0, 0)]);

    x[IX(0, N - 1)] = 0.33f * (x[IX(1, N - 1)]
        + x[IX(0, N - 2)]
        + x[IX(0, N - 1)]);

    x[IX(N - 1, 0)] = 0.33f * (x[IX(N - 2, 0)]
        + x[IX(N - 1, 1)]
        + x[IX(N - 1, 0)]);

    x[IX(N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1)]
        + x[IX(N - 1, N - 2)]
        + x[IX(N - 1, N - 1)]);
};

static void lin_solve(int b, std::vector<float> x, std::vector<float> x0, float a, float c) {
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

static void diffuse(int b, std::vector<float> x, std::vector<float> x0, float diff, float dt) {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 4 * a);
};

static void project(std::vector<float> velocX, std::vector<float> velocY, std::vector<float> p, std::vector<float> div) {
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
    lin_solve(0, p, div, 1, 6);

    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] = -0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] = -0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
        }
    }

    set_bnd(1, velocX);
    set_bnd(2, velocY);
};

static void advect(int b, std::vector<float> d, std::vector<float> d0, std::vector<float> velocX, std::vector<float> velocY, float dt) {
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
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
            i0 = floorf(x);
            i1 = i0 + 1.0f;
            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
            j0 = floorf(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;

            d[IX(i, j)] =
                s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
        }
    }

    set_bnd(b, d);
};

void step(fluidGridType* grid) {
    float visc = grid->visc;
    float diff = grid->diff;
    float dt = grid->dt;

    std::vector<float> Vx = grid->Vx;
    std::vector<float> Vy = grid->Vy;
    std::vector<float> Vx0 = grid->Vx0;
    std::vector<float> Vy0 = grid->Vy0;

    std::vector<float> s = grid->s;
    std::vector<float> density = grid->density;

    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);

    project(Vx0, Vy0, Vx, Vy);

    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);

    project(Vx, Vy, Vx0, Vy0);

    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
};


void draw(fluidGridType* grid, unsigned int texture) {
    //step(grid);

    std::vector<unsigned char> textureData(N * N * 4);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            textureData[4 * IX(i, j)] = 255.0;  // R
            textureData[4 * IX(i, j) + 1] = 255.0;  // G
            textureData[4 * IX(i, j) + 2] = 255.0;  // B
            textureData[4 * IX(i, j) + 3] = grid->density[IX(i, j)];
            //grid->density[IX(i, j)] += 5.0;
            //std::cout << "grid density " << grid->density[IX(i, j)] << std::endl;
            //textureData[IX(i, j) + 3] = (grid->density[IX(i, j)]);
            //if (grid->density[IX(i, j)] != 0) {
            //    std::cout << "val " << (grid->density[IX(i, j)]) << std::endl;
            //   // std::cout << "Index " << IX(i, j) << std::endl;
            //   // std::cout << "grid density " << grid->density[IX(i, j)] << std::endl;
            //}
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


    float vertices1[] = {
        -0.5f, -0.5f, 0.0f,
        0.5f, -0.5f, 0.0f,
        -0.0f, 0.5f, 0.0f,
    };

    float vertices2[] = {
        -0.9f, -0.9f, 0.0f,
        0.4f, -0.2f, 0.0f,
        -0.0f, 0.5f, 0.0f,
    };

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

        // Set the uniform variables
        //glUniform2f(glGetUniformLocation(shaderProgram, "uScreenSize"), 640, 480);
        //glUniform2f(glGetUniformLocation(shaderProgram, "uPosition"), rectX, rectY);
        //glUniform2f(glGetUniformLocation(shaderProgram, "uSize"), rectWidth, rectHeight);

        glBindVertexArray(VAO);

        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        //glDrawArrays(GL_TRIANGLES, 0, 3);
        //glDrawArrays(GL_TRIANGLES, 3, 6);

        //draw(&fluidGrid, shaderProgram);

        //glUseProgram(program2);

        //glBindVertexArray(VAO[1]);
        //glDrawArrays(GL_TRIANGLES, 0, 3);

        //glBegin(GL_TRIANGLES);
        //glVertex2f(-0.5f, -0.5f);
        //glVertex2f(0.0f, 0.5f);
        //glVertex2f(0.5f, -0.5f);
        //glEnd();

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

        //float textureX = (rectX / SCREEN_WIDTH) * SCREEN_WIDTH;
        //float textureY = (rectY / SCREEN_HEIGHT) * SCREEN_HEIGHT;

        //// Update rectX and rectY to be within the bounds of the texture
        //rectX = textureX;
        //rectY = textureY;
        fluidGridType* fluidGrid = static_cast<fluidGridType*>(glfwGetWindowUserPointer(window));
        std::cout << "rect x adn y " << rectX << " " << rectY << " " << std::endl;
        addDensity(fluidGrid, rectX, rectY, 100.0);
        //addVelocity(fluidGrid, rectX, rectY, 5.0, 5.0);
    }
}