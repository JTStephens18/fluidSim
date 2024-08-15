#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);

static const char* vertexShaderSource = "#version 460 core\n"
"layout (location = 0) in vec3 aPos; \n"
"void main()\n"
"{\n"
"gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
"}\0";

static const char* fragmentShaderSource = "#version 460 core\n"
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

int N = 64;
int iter = 4;

int IX(int x, int y) {
    return ((x)+N * (y));
};

struct fluidGrid {

    int N;

    std::vector<float> Vx;
    std::vector<float> Vy;

    std::vector<float> Vx0;
    std::vector<float> Vy0;

    std::vector<float> density;
    std::vector<float> s;
};

fluidGrid fluidGridCreate() {
    fluidGrid newGrid;
    newGrid.N = N;
    int gridSize = newGrid.N * newGrid.N;
    newGrid.Vx.assign(gridSize, 0.0);
    newGrid.Vy.assign(gridSize, 0.0);
    newGrid.Vx0.assign(gridSize, 0.0);
    newGrid.Vy0.assign(gridSize, 0.0);
    newGrid.density.assign(gridSize, 0.0);
    newGrid.s.assign(gridSize, 0.0);
    return newGrid;
};

void addDensity(fluidGrid* grid, int x, int y, float amount) {
    grid->density[IX(x, y)] += amount;
};

void addVelocity(fluidGrid* grid, int x, int y, float amountX, float amountY) {
    int index = IX(x, y);
    grid->Vx[index] += amountX;
    grid->Vy[index] += amountY;
};

static void lin_solve

static void diffuse(int b, float* x, float* x0, float diff, float dt) {
    float a = dt * diff * (N - 2) * (N - 2);
    //lin_solve(b, x, x0, a, 1 + 4 * a);
};


int main(void)

{
    fluidGrid fluidGrid = fluidGridCreate();



    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    /*
    glfwCreateWindow takes in width, height, name
        returns a GLFWwindow object
    */
    window = glfwCreateWindow(640, 480, "Hello World Triangle", NULL, NULL);
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
    /* We register callback functions after creating the window, but before the render loop */
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    unsigned int vertexShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    unsigned int fragmentShader;
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    unsigned int fragmentShader2;
    fragmentShader2 = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader2, 1, &fragmentShaderSource2, NULL);
    glCompileShader(fragmentShader2);

    int success;
    char infoLog[512];
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    //glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" <<
            infoLog << std::endl;
    }

    unsigned int shaderProgram, program2;
    shaderProgram = glCreateProgram();
    program2 = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glAttachShader(program2, vertexShader);
    glAttachShader(program2, fragmentShader2);
    glLinkProgram(program2);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glDeleteShader(fragmentShader2);

    float vertices[] = {
        -0.5f, -0.5f, 0.0f,
        0.5f, -0.5f, 0.0f,
        -0.0f, 0.5f, 0.0f,
    };

    float vertices2[] = {
        -0.9f, -0.9f, 0.0f,
        0.4f, -0.2f, 0.0f,
        -0.0f, 0.5f, 0.0f,
    };

    unsigned int VBO[2], VAO[2];
    /* A vertex array object (VAO) can be bound like a vertex buffer object. It has the advantage of storting vertex attribute pointers (such as vertex shaders) so we only have to make the calls once */
    glGenVertexArrays(2, VAO);
    glGenBuffers(2, VBO);

    glBindVertexArray(VAO[0]);
    glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(VAO[1]);
    glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices2), vertices2, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //glBindBuffer(GL_ARRAY_BUFFER, 0);


    // Unbind VAO
    //glBindVertexArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {

        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);

        glBindVertexArray(VAO[0]);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        //glDrawArrays(GL_TRIANGLES, 3, 6);

        glUseProgram(program2);

        glBindVertexArray(VAO[1]);
        glDrawArrays(GL_TRIANGLES, 0, 3);

        //glBegin(GL_TRIANGLES);
        //glVertex2f(-0.5f, -0.5f);
        //glVertex2f(0.0f, 0.5f);
        //glVertex2f(0.5f, -0.5f);
        glEnd();

        processInput(window);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events (such as inputs) */
        glfwPollEvents();
    }
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