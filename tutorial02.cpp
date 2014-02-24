// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <ctime>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <glfw3.h>
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace glm;

// Include AntTweakBar
#include <AntTweakBar.h>

#include <common/shader.hpp>
#include <math.h>
using namespace glm;
using namespace std; 
#include <iostream>

vec3 gPosition1(-1.5f, 0.0f, 0.0f);

float euler(float stepSize, float lastValue, float yprim) {
    return lastValue + stepSize*yprim;
}

int main( int argc, char **argv )
{
	// Initialise GLFW
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow( 1024, 768, "Tutorial 02 - Red triangle", NULL, NULL);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// Dark blue background
	glClearColor(0.8f, 0.9f, 1.0f, 0.0f);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders( "SimpleVertexShader.vertexshader", "SimpleFragmentShader.fragmentshader" );

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");

	// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	glm::mat4 Projection = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);
	// Camera matrix
	glm::mat4 View       = glm::lookAt(
								glm::vec3(0,0,-3), // Camera is at (4,3,-3), in World Space
								glm::vec3(0,0,0), // and looks at the origin
								glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
						   );
	// Model matrix : an identity matrix (model will be at the origin)
	glm::mat4 Model      = glm::mat4(1.0f);
	// Our ModelViewProjection : multiplication of our 3 matrices
	glm::mat4 MVP        = Projection * View * Model; // Remember, matrix multiplication is the other way around
	
	const int n_ballz = 7;
	GLfloat ballz[190*n_ballz] = {0.0f};

	GLfloat colorz[189*n_ballz] = {0.0f};
	float ballcolors[3*n_ballz]= {0.0f};
	GLuint colorBallBuffer; 

	
	srand(time(NULL));

	for(int c=0; c<3*n_ballz; c+=3)
	{
		ballcolors[c] = (double)(rand() % 100)/100;
		ballcolors[c+1] = (double)(rand() % 100)/100;
		ballcolors[c+2] = (double)(rand() % 100)/100;
	}


	for(int i=0; i<189*n_ballz; i= i+3)
	{
		int j = i/189;
		colorz[i] = ballcolors[j*3];
		colorz[i+1] = ballcolors[j*3+1];
		colorz[i+2] = ballcolors[j*3+2];
	}


	//wallz
	float wlength = 1.0f;
	float wheight = 0.7f;
	float wwidth = 0.1f;

	static const GLfloat walls[] = {
	  //West wall
	  -wlength, wheight, 0.0f,
	  -wlength + wwidth, wheight, 0.0f,
	  -wlength + wwidth, -wheight, 0.0f,

	  -wlength + wwidth, -wheight, 0.0f,
	  -wlength, -wheight, 0.0f,
	  -wlength, wheight, 0.0f,

	  //South wall
	  -wlength, -wheight, 0.0f,
	  wlength, -wheight, 0.0f,
	  wlength, -wheight - wwidth, 0.0f,

	  wlength, -wheight - wwidth, 0.0f,
	  -wlength, -wheight - wwidth, 0.0f,
	  -wlength, -wheight, 0.0f,

	  //East wall
	  wlength - wwidth, wheight, 0.0f,
	  wlength, wheight, 0.0f,
	  wlength, -wheight, 0.0f,

	  wlength, -wheight, 0.0f,
	  wlength - wwidth, -wheight, 0.0f,
	  wlength - wwidth, wheight, 0.0f,
	};
	
	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(walls), walls, GL_STATIC_DRAW);

	GLuint ropebuffer; 
	GLuint ballzbuffer;


	//Const
	float pi = 3.14159f, stepSize = 0.005f, g = 9.82f;
	bool isPressed = false;
	bool firstCheck = false;

	//Initialvalues
	float theta[n_ballz];
	float velocity[n_ballz];
	float acceleration[n_ballz];
	float ropeLength[n_ballz];

	float xAcc[n_ballz] = {0.0f};
	float yAcc[n_ballz] = {0.0f};
	float xVel[n_ballz] = {0.0f};
	float yVel[n_ballz] = {0.0f};
	float xAirres[n_ballz] = {0.0f};
	float yAirres[n_ballz] = {0.0f};
	float airres[n_ballz] = {0.0f};
	float xPosition[n_ballz] = {0.0f};
	float yPosition[n_ballz] = {0.0f};
	float prexPos[n_ballz] = {0.0f};
	float preyPos[n_ballz] = {0.0f};
	float prexVel[n_ballz] = {0.0f};
	float preyVel[n_ballz] = {0.0f};

	float radius[n_ballz] = {0.0f};
	float volume[n_ballz] = {0.0f};
	float area[n_ballz] = {0.0f};
	float density[n_ballz] = {0.0f}; 
	float mass[n_ballz] = {0.0f};
	float airconstant[n_ballz] = {0.0f};
	//Balls
	for(int l = 0; l<n_ballz; l++){
		radius[l] = 0.03f;
		volume[l] = ((4.0f*pi*pow(radius[l],3.0f))/3.0f);
		area[l] = 4.0f*pi*pow(radius[l],2.0f);
		density[l] = 11340.0f; 
		mass[l] = density[l]*volume[l];
		airconstant[l] = 0.47f;
		theta[l] = pi/2;
		velocity[l] = 0.0f;
		acceleration[l] = 0.0f;
		ropeLength[l] = 0.2f;
	}
	/* 
	//Om man vill välja lite values själv :)
	for(int l = 0; l<n_ballz; l++){

		float degree;

		cout << "Boll: "  << (l+1) << endl;
		cout << "Radie: ";
		cin >> radius[l];
		cout << "Startvinkel i grader: ";
		cin >> degree;
		theta[l] = degree*3.14/180;
		cout << "Vinkelhastighet: ";
		cin >> velocity[l];
		cout << "Length (0.05 - 0.2): ";
		cin >> ropeLength[l];
		volume[l] = ((4.0f*pi*pow(radius[l],3.0f))/3.0f);
		area[l] = 4.0f*pi*pow(radius[l],2.0f);
		density[l] = 11340.0f; 
		mass[l] = density[l]*volume[l];
		airconstant[l] = 0.47f;
		acceleration[l] = 0.0f;
	}*/

	// Calc ropes
	const float distBetw = 0.1f;
	GLfloat ropez[6*n_ballz] = {0.0f};

	//If uneven n_ballz
	if(n_ballz==1){
		ropez[0] = 0;
		ropez[1] = ropeLength[0];
	}
	
	else if(n_ballz%2 == 1){
		for(int k = 0; k<n_ballz; k++) {
			if (k==0){
				ropez[k]=0;
				ropez[k+1] = ropeLength[k];
			}
			//Negative
			else if(k%2 == 1) {
				ropez[k*6] = -(k+1)*distBetw;
				ropez[k*6+1] = ropeLength[k]; 
			}
			//Positive
			else {
				ropez[k*6] = -ropez[(k-1)*6];
				ropez[k*6+1] = ropeLength[k]; 
			}
		}
	}
	//If even n_ballz
	else{
		for(int k = 0; k<n_ballz; k++) {
			//Negative
			if(k%2 != 1) {
				ropez[k*6] = (k+1)*distBetw;
				ropez[k*6+1] = ropeLength[k]; 
			}
			//Positive
			else {
				ropez[k*6] = -ropez[(k-1)*6];
				ropez[k*6+1] = ropeLength[k]; 
			}
		}
	}

	// For speed computation
	double lastTime = glfwGetTime();
	double lastFrameTime = lastTime;
	float timeBetween = 0.0f;

	do{
		// Measure speed
		double currentTime = glfwGetTime();
		float deltaTime = (float)(currentTime - lastFrameTime); 
		timeBetween += deltaTime;
		lastFrameTime = currentTime;

		if (timeBetween >= 1/30 ){ // calc each 1/60 sec
			timeBetween = 0;

			// Clear the screen
			glClear( GL_COLOR_BUFFER_BIT );

			// Use our shader
			glUseProgram(programID);

			// Send our transformation to the currently bound shader, 
			// in the "MVP" uniform
			glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

			//If space is pressed; throw ball!
			if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && isPressed == false) {
				isPressed = true;
				for(int r=0; r<sizeof(ropez)/sizeof(*ropez); r++)
					ropez[r] = 0.0f;
				
				glGenBuffers(1, &ropebuffer); 
				glBindBuffer(GL_ARRAY_BUFFER, ropebuffer); 
				glBufferData(GL_ARRAY_BUFFER, sizeof(ropez), ropez, GL_STATIC_DRAW);

			}
			//Pendulum
			if (isPressed == false) {
				for(int index = 0; index<n_ballz; index++) 
				{
					airres[index] = (0.5f*pow(velocity[index],2.0f)*area[index]*airconstant[index])/mass[index];
					

					if (velocity[index] < 0) 
						airres[index] = airres[index]*-1.0f;
				
					acceleration[index] = (-g / ropeLength[index])*sin(theta[index]);
					velocity[index] = euler(stepSize, velocity[index], acceleration[index]) - airres[index];
					theta[index] = euler(stepSize, theta[index], velocity[index]);
					xPosition[index] = ropeLength[index]*sin(theta[index]);
					yPosition[index] = ropeLength[index]*(1 - cos(theta[index]));

					xPosition[index] += ropez[index*6];

					ropez[index*6 + 3] = xPosition[index];
					ropez[index*6 + 4] = yPosition[index];

					glGenBuffers(1, &ropebuffer); 
					glBindBuffer(GL_ARRAY_BUFFER, ropebuffer); 
					glBufferData(GL_ARRAY_BUFFER, sizeof(ropez), ropez, GL_STATIC_DRAW);
				}
			}
    
			//Projectile
			else {
				if (firstCheck == false){
					for(int index = 0; index<n_ballz; index++) {
						xAcc[index] = 0.0f; 
						yAcc[index] = acceleration[index] * sin(theta[index]) - g; 
						xVel[index] = velocity[index] * cos(theta[index]); 
						yVel[index] = velocity[index] * sin(theta[index]);
						ropez[index*6 +1] = 0.0f;
						ropez[index*6 +3] = 0.0f;
						ropez[index*6 +4] = 0.0f;
						glGenBuffers(1, &ropebuffer); 
						glBindBuffer(GL_ARRAY_BUFFER, ropebuffer); 
						glBufferData(GL_ARRAY_BUFFER, sizeof(ropez), ropez, GL_STATIC_DRAW);
					}
					firstCheck = true;
				}
				for(int index = 0; index < n_ballz; index++) {
					xAirres[index] = (0.5f*pow(xVel[index],2)*area[index]*airconstant[index])/mass[index];
					yAirres[index] = (0.5f*pow(yVel[index],2)*area[index]*airconstant[index])/mass[index];
        
					if (xVel[index] < 0) 
						xAirres[index] = -1*xAirres[index];
				

					if (yVel[index] < 0) 
						yAirres[index] = -1*yAirres[index];
				
					xVel[index] = euler(stepSize, xVel[index], xAcc[index]) - xAirres[index];
					yVel[index] = euler(stepSize, yVel[index], yAcc[index]) - yAirres[index];
					xPosition[index] = euler(stepSize, xPosition[index], xVel[index]);
					yPosition[index] = euler(stepSize, yPosition[index], yVel[index]);

					//Hit ground (add radius)
					if (yPosition[index] - radius[index] < -wheight){
						yVel[index] = -0.5*yVel[index];
						yPosition[index] = -wheight + radius[index];
						xVel[index] = xVel[index]*0.5f;
					}
					//Hit the west or east wall (add radius)
					if (xPosition[index] - radius[index] < -wlength + wwidth) {
						xVel[index] = -0.5*xVel[index];
						xPosition[index] = -wlength + wwidth + radius[index];
					}
					if (xPosition[index] + radius[index] > wlength - wwidth) {
						xVel[index] = -0.5*xVel[index];
						xPosition[index] = wlength - wwidth - radius[index];
					}
					ballz[index*3] = xPosition[index]; 
					ballz[index*3+1] = yPosition[index];
				}

				
				//Collision between two bouncing objects      
				bool collision;
				do
				{
					collision = false;
					for(int one = 0; one < n_ballz; one++)
					{
						for(int two = one+1; two < n_ballz; two++)
						{
							if(sqrtf(pow(xPosition[one]-xPosition[two], 2) + pow(yPosition[one]-yPosition[two], 2)) < radius[one] + radius[two] )
							{
								cout << "collision!! "<< endl;
								collision = true;
								xVel[one] = prexVel[one];
								xVel[two] = prexVel[two]; 

								yVel[one] = preyVel[one];
								yVel[two] = preyVel[two];

							
								xPosition[one] = prexPos[one];
								xPosition[two] = prexPos[two];

								yPosition[one] = preyPos[one];
								yPosition[two] = preyPos[two];

								float angleOne = atan2(yVel[one], xVel[one]);
								float angleTwo = atan2(yVel[two], xVel[two]);
							
								float move = 0.001;
								while(sqrtf(pow(xPosition[one]-xPosition[two], 2) + pow(yPosition[one]-yPosition[two], 2)) < radius[one] + radius[two])
								{								
									float tempXPos[] = {xPosition[one], xPosition[two]}; 
									float tempYPos[] = {yPosition[one], yPosition[two]}; 

									xPosition[one] = xPosition[one] - (move * cos(angleOne));
									yPosition[one] = yPosition[one] - (move * sin(angleOne));

									xPosition[two] = xPosition[two] - (move * cos(angleTwo));
									yPosition[two] = yPosition[two] - (move * sin(angleTwo));

									//Ground
									if (yPosition[one] - radius[one] < -wheight)
										yPosition[one] = -wheight + radius[one];
									if (yPosition[two] - radius[two] < -wheight)
										yPosition[two] = -wheight + radius[two];
								
									//West Wall
									if (xPosition[one] - radius[one] < -wlength + wwidth)
										xPosition[one] = -wlength + wwidth + radius[one];
									if (xPosition[two] - radius[two] < -wlength + wwidth)
										xPosition[two] = -wlength + wwidth + radius[two];
									//West Wall
									if (xPosition[one] + radius[one] > wlength - wwidth)
										xPosition[one] = wlength - wwidth -radius[one];
									if (xPosition[two] + radius[two] > wlength - wwidth)
										xPosition[two] = wlength - wwidth -radius[two];

									if(sqrtf(pow(xPosition[one]-xPosition[two], 2) + pow(yPosition[one]-yPosition[two], 2)) > radius[one] + radius[two])
									{
										xPosition[one] = tempXPos[0];
										yPosition[one] = tempYPos[0];

										xPosition[two] = tempXPos[1];
										yPosition[two] = tempYPos[1];
										break;
									}

								}

								vec2 posOne = vec2(xPosition[one], yPosition[one]);
								vec2 posTwo = vec2(xPosition[two], yPosition[two]);
								vec2 xPos = posTwo - posOne;

								vec2 x = normalize(xPos);
								vec2 v1 = vec2(xVel[one], yVel[one]);
								float x1 = dot(x, v1);
								vec2 v1x = x*x1; 
								vec2 v1y = v1 - v1x; 

								x = -x; 
								vec2 v2 = vec2(xVel[two], yVel[two]);
								float x2 = dot(x, v2); 
								vec2 v2x = x*x2; 
								vec2 v2y = v2 -v2x; 

								float totMass = mass[one]+mass[two];

								vec2 newV1 = (v1x*((mass[one]-mass[two]) /totMass)) + (v2x*(mass[two]/totMass)) + v1y;
								vec2 newV2 = (v1x*((mass[one] /totMass)) + (v2x*(mass[two]-mass[one])/totMass)) + v2y;

								xVel[one] = newV1[0]*0.8f;
								xVel[two] = newV1[1]*0.8f;

								yVel[one] = newV2[0]*0.8f;
								yVel[two] = newV2[1]*0.8f;
							}
						}
					}
				}while(collision == true);
			}

			//draw the two ballz
			for(int index = 0; index < n_ballz; index++){
				int i = 189*index;
				float ang = 0.0f;
				float step = 0.1f; //same decimals for even size of g_vertex_buffer_data[int]
				float pi2 = 3.1f;

				//cout << index << " x: " << xVel[index]  << " y: " << yPosition[index] << endl;
      
				do{
					ballz[i] = radius[index]*cos(ang) + xPosition[index];
					i++;
					ballz[i] = radius[index]*sin(ang) + yPosition[index];
					i += 2;
					ang += step;
				} while (ang <= (2 * pi2));
    
				ballz[i] = radius[index] + xPosition[index];
				ballz[i+1] = yPosition[index];
			} 
			for(int g = 0; g < n_ballz; g++)
			{
				prexPos[g] = xPosition[g];
				preyPos[g] = yPosition[g];
				prexVel[g] = xVel[g];
				preyVel[g] = yVel[g];
			}

			glGenBuffers(1, &ballzbuffer); 
			glBindBuffer(GL_ARRAY_BUFFER, ballzbuffer); 
			glBufferData(GL_ARRAY_BUFFER, sizeof(ballz), ballz, GL_STATIC_DRAW);

			glGenBuffers(1, &colorBallBuffer); 
			glBindBuffer(GL_ARRAY_BUFFER, colorBallBuffer); 
			glBufferData(GL_ARRAY_BUFFER, sizeof(colorz), colorz, GL_STATIC_DRAW);


			int width, height;
			glfwGetFramebufferSize(window, &width, &height);
			glViewport(0, 0, width, height);
			glClear(GL_COLOR_BUFFER_BIT);

			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
			glVertexAttribPointer(
									0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
									3,                  // size
									GL_FLOAT,           // type
									GL_FALSE,           // normalized?
									0,                  // stride
									(void*)0            // array buffer offset
									);

			// Draw the triangle !
			glDrawArrays(GL_TRIANGLES, 0, 3*6); // 3 indices starting at 0 -> 1 triangle

			//for the rope
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, ropebuffer);
			glVertexAttribPointer(
									0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
									3,                  // size
									GL_FLOAT,           // type
									GL_FALSE,           // normalized?
									0,                  // stride
									(void*)0            // array buffer offset
									);
    

			//glVertexPointer(2, GL_FLOAT, 0, ropez);
			glDrawArrays(GL_LINES, 0, 2*n_ballz);
			glEnableVertexAttribArray(1); 
    
			//for them baallz
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, ballzbuffer);
			glVertexAttribPointer(
									0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
									3,                  // size
									GL_FLOAT,           // type
									GL_FALSE,           // normalized?
									0,                  // stride
									(void*)0            // array buffer offset
									);

			glBindBuffer(GL_ARRAY_BUFFER, colorBallBuffer); 
			glVertexAttribPointer(
									1,                  // attribute 1. No particular reason for 0, but must match the layout in the shader.
									3,                  // size
									GL_FLOAT,           // type
									GL_FALSE,           // normalized?
									0,                  // stride
									(void*)0            // array buffer offset
									);

			// Draw the triangle !
			for(int w = 0; w < n_ballz; w++)
				glDrawArrays(GL_TRIANGLE_FAN, w*63, 63); // 3 indices starting at 0 -> 1 triangle
	

			//Disable
			glDisableVertexAttribArray(1); 
			glDisableVertexAttribArray(0);
			
			// Swap buffers
			glfwSwapBuffers(window);
			glfwPollEvents();
		}
	} // Check if the ESC key was pressed or the window was closed
	while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 );

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	// Cleanup VBO
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteVertexArrays(1, &VertexArrayID);

	return 0;
}

