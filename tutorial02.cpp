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

// Include AntTweakBar
#include <AntTweakBar.h>

#include <common/shader.hpp>
#include <math.h>
using namespace glm;
using namespace std; 
#include <iostream>

//numerical integration method: euler
float euler(float stepSize, float lastValue, float yprim) {
    return lastValue + stepSize*yprim;
}
int main( int argc, char **argv )
{
	// Initialise GLFW, with errormsg
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

	//background colour
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
								glm::vec3(0,0,-5), // Camera is at (4,3,-3), in World Space
								glm::vec3(0,0,0), // and looks at the origin
								glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
						   );
	// Model matrix : an identity matrix (model will be at the origin)
	glm::mat4 Model      = glm::mat4(1.0f);
	// Our ModelViewProjection : multiplication of our 3 matrices
	glm::mat4 MVP        = Projection * View * Model; // Remember, matrix multiplication is the other way around
	
	// How many balls we want in the simulation
	const int n_ballz = 15;
	GLfloat ballz[190*n_ballz] = {0.0f};
	
	// Initialize color of the balls
	GLfloat colorz[189*n_ballz] = {0.0f};
	float ballcolors[3*n_ballz]= {0.0f}; //each ball get an own color (r,g,b)
	GLuint colorBallBuffer; 
	srand(time(NULL));

	// Set every ball do a random color
	for(int c=0; c<3*n_ballz; c+=3)
	{
		ballcolors[c] = (double)(rand() % 100)/100;
		ballcolors[c+1] = (double)(rand() % 100)/100;
		ballcolors[c+2] = (double)(rand() % 100)/100;
	}
	//apply the color each ball
	for(int i=0; i<189*n_ballz; i= i+3)
	{
		int j = i/189;
		colorz[i] = ballcolors[j*3];
		colorz[i+1] = ballcolors[j*3+1];
		colorz[i+2] = ballcolors[j*3+2];
	}


	//Set how big the walls should be
	float wlength = 2.0f;
	float wheight = 0.7f;
	float wwidth = 0.1f; //of the walls and floor

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
	//This is the wallsbuffer
	GLuint vertexbuffer;
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(walls), walls, GL_STATIC_DRAW);

	//Initialize the rope- and ballbuffer;
	GLuint ropebuffer; 
	GLuint ballzbuffer;

	//Initialize som constants
	float pi = 3.14159f, stepSize = 0.005f, g = 9.82f;
	bool isPressed = false;
	bool firstCheck = false;

	//Allocate arrays for the objects' properties
	float theta[n_ballz];
	float velocity[n_ballz];
	float acceleration[n_ballz];
	float ropeLength[n_ballz];

	//initialize them
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

	//Initialize for non-zero values
	for(int l = 0; l<n_ballz; l++){
		radius[l] = 0.05f;
		volume[l] = ((4.0f*pi*pow(radius[l],3.0f))/3.0f);
		area[l] = 4.0f*pi*pow(radius[l],2.0f);
		density[l] = 11340.0f; 
		mass[l] = density[l]*volume[l];
		airconstant[l] = 0.47f;
		theta[l] = pi/2;
		velocity[l] = 30.0f;
		acceleration[l] = 0.0f;
		ropeLength[l] = 0.2f;
	}

	/*
	//Om man vill välja lite values själv :)
	for(int l = 0; l<n_ballz; l++){

		float degree;

		cout << "Boll: "  << (l+1) << endl;
		cout << "Radie (0.01 - 0.05): ";
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

	/***********************************************************/
	// Calc how the balls(ropes) should appear in the simulation
	/***********************************************************/
	const float distBetw = 0.1f; //between two penduluma
	GLfloat ropez[6*n_ballz] = {0.0f};

	//to draw the pendulums within nice intervals and distance to eachother and the walls
	// If it is only one ball
	if(n_ballz==1){
		ropez[0] = 0;
		ropez[1] = ropeLength[0];
	}
	//If uneven number of balls
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
	//If even number of balls
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

	//the render loop!
	do{

		// Clear the screen
		glClear( GL_COLOR_BUFFER_BIT );

		// Use our shader
		glUseProgram(programID);

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

		//If space abar is pressed; throw ball!
		if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS && isPressed == false) {
			isPressed = true;
			//the ropes should be erased and their buffer binded
			for(int r=0; r<sizeof(ropez)/sizeof(*ropez); r++)
				ropez[r] = 0.0f;
			
			glGenBuffers(1, &ropebuffer); 
			glBindBuffer(GL_ARRAY_BUFFER, ropebuffer); 
			glBufferData(GL_ARRAY_BUFFER, sizeof(ropez), ropez, GL_STATIC_DRAW);

		}
		//if the space bar is not pressed - the balls should be drawn as pendulums
		if (isPressed == false) {
			for(int index = 0; index<n_ballz; index++) 
			{
				airres[index] = (0.5f*pow(velocity[index],2.0f)*area[index]*airconstant[index])/mass[index];
					
				//the air resistant is always opposite to the velocity
				if (velocity[index] < 0) 
					airres[index] = airres[index]*-1.0f;
				
				/* Using Euler
				acceleration[index] = (-g / ropeLength[index])*sin(theta[index]);
				velocity[index] = euler(stepSize, velocity[index], acceleration[index]) - airres[index];
				theta[index] = euler(stepSize, theta[index], velocity[index]);*/

				// Using Runge-kutta4
				float k1theta = velocity[index];
				float k1w = -g/ropeLength[index]*sin(theta[index]);

				float k2theta = velocity[index]+stepSize/2*k1w;
				float k2w = (-g/ropeLength[index])*sin(theta[index]+k1theta*stepSize/2);

				float k3theta = velocity[index]+stepSize/2*k2w;
				float k3w = (-g/ropeLength[index])*sin(theta[index]+k2theta*stepSize/2);

				float k4theta = velocity[index]+stepSize*k2w;
				float k4w = (-g/ropeLength[index])*sin(theta[index]+k2theta*stepSize);

				//get acceleration, velocity and position for all balls
				acceleration[index] = -(g/ropeLength[index]*sin(theta[index]));
				velocity[index] = velocity[index]+(stepSize*(k1w + 2*k2w + 2*k3w + k4w)/6)-airres[index];
				theta[index] = theta[index]+(stepSize*(k1theta + 2*k2theta + 2*k3theta + k4theta)/6);

				xPosition[index] = ropeLength[index]*sin(theta[index]);
				yPosition[index] = ropeLength[index]*(1 - cos(theta[index]));

				xPosition[index] += ropez[index*6];
				//change the position of the rope according to its balls xposition
				ropez[index*6 + 3] = xPosition[index];
				ropez[index*6 + 4] = yPosition[index];

				//bind the buffers
				glGenBuffers(1, &ropebuffer); 
				glBindBuffer(GL_ARRAY_BUFFER, ropebuffer); 
				glBufferData(GL_ARRAY_BUFFER, sizeof(ropez), ropez, GL_STATIC_DRAW);
			}
		}
    
		//If the space bar is pressed: Projectile movement
		else {
			if (firstCheck == false){
				for(int index = 0; index<n_ballz; index++) {
					xAcc[index] = 0.0f; 
					yAcc[index] = acceleration[index] * sin(theta[index]) - g; 
					xVel[index] = velocity[index] * cos(theta[index]);
					yVel[index] = velocity[index] * sin(theta[index]);

					//in order to erase the ropes -> zero values. 
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
				/*using eluer
				xVel[index] = euler(stepSize, xVel[index], xAcc[index]) - xAirres[index];
				yVel[index] = euler(stepSize, yVel[index], yAcc[index]) - yAirres[index];
				xPosition[index] = euler(stepSize, xPosition[index], xVel[index]);
				yPosition[index] = euler(stepSize, yPosition[index], yVel[index]);*/

				//Using RungeKutta4
				float k1xPos = xVel[index];
				float k1yPos = yVel[index];
				float k1xVel = 0;
				float k1yVel = -g;
					
				float k2xPos = xVel[index]+stepSize/2*k1xVel;
				float k2yPos = yVel[index]+stepSize/2*k1yVel;
				float k2yVel = -g + k1yPos*stepSize/2;
				float k2xVel = k1xPos*stepSize/2;

				float k3xPos = xVel[index]+stepSize/2*k2xVel;
				float k3yPos = yVel[index]+stepSize/2*k2yVel;
				float k3yVel = -g + k2yPos*stepSize/2;
				float k3xVel = k2xPos*stepSize/2;

				float k4xPos = xVel[index]+stepSize*k3xVel;
				float k4yPos = yVel[index]+stepSize*k3yVel;
				float k4yVel = -g + k3yPos*stepSize;
				float k4xVel = k3xPos*stepSize;

				yPosition[index] = yPosition[index] + (stepSize*(k1yPos + 2*k2yPos + 2*k3yPos + k4yPos)/6);
				xPosition[index] = xPosition[index] + (stepSize*(k1xPos + 2*k2xPos + 2*k3xPos + k4xPos)/6);
				yVel[index] = yVel[index] + (stepSize*(k1yVel + 2*k2yVel + 2*k3yVel + k4yVel)/6);
				xVel[index] = xVel[index] + (stepSize*(k1xVel + 2*k2xVel + 2*k3xVel + k4xVel)/6);

				//collision detection!
				//1. check the distance to the ground 
				if (yPosition[index] - radius[index] < -wheight){ 	//if the edge of the circle is below the ground, 
					yVel[index] = -0.5*yVel[index];					//change the direction of the velocity and damp with x0.5
					xVel[index] = 0.5f*xVel[index];
					yPosition[index] = -wheight + radius[index];	//set the x-pos of the ball to lie on top of the floor
					
				}
				//2. Hit the west wall
				if (xPosition[index] - radius[index] < -wlength + wwidth) {	
					xVel[index] = -0.5*xVel[index];							//change direction of the x-velocity and damp with a factor of 0.5
					xPosition[index] = -wlength + wwidth + radius[index];	//set the x-position of the ball to lie adjacent to the west wall
				}
				//3. Hit the east wall
				if (xPosition[index] + radius[index] > wlength - wwidth) {	
					xVel[index] = -0.5*xVel[index];							//change direction of the x-velocity and damp with a factor of 0.5
					xPosition[index] = wlength - wwidth - radius[index];	//set the x-position of the ball to lie adjacent to the east wall
				}
				//set the values for the ball!
				ballz[index*3] = xPosition[index]; 
				ballz[index*3+1] = yPosition[index];
			}

			//Collision between two bouncing objects      
			bool collision;
			do
			{
				//Keep going 
				collision = false;
				for(int one = 0; one < n_ballz; one++) //check collision between two balls. 
				{
					for(int two = one+1; two < n_ballz; two++)
					{
						//if the distance(hypotenuse) between the balls are smaller than the sum of their radius
						//they have stuck together -> fix so that they are two separate objects
						if(sqrtf(pow(xPosition[one]-xPosition[two], 2) + pow(yPosition[one]-yPosition[two], 2)) < radius[one] + radius[two] )
						{	
							//save the previous values of the balls
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

							//calculate the normalplane
							vec2 posOne = vec2(xPosition[one], yPosition[one]);
							vec2 posTwo = vec2(xPosition[two], yPosition[two]);
							vec2 norm = posOne - posTwo;
							vec2 normalPlane = normalize(norm);
							//move the ball from the previous values: along the normal plane until they are adjacents
							while(sqrtf(pow(xPosition[one]-xPosition[two], 2) + pow(yPosition[one]-yPosition[two], 2)) < radius[one] + radius[two])
							{			
									
								float tempXPos[] = {xPosition[one], xPosition[two]}; 
								float tempYPos[] = {yPosition[one], yPosition[two]}; 

								//move towards each other: in the right direction
								if(xPosition[one] < xPosition[two]){
									posOne -= move*normalPlane;
									posTwo += move*normalPlane;
								}
								else{
									posOne += move*normalPlane;
									posTwo -= move*normalPlane;
								}

								//Check for collision with the ground
								if (yPosition[one] - radius[one] < -wheight)
									yPosition[one] = -wheight + radius[one];
								if (yPosition[two] - radius[two] < -wheight)
									yPosition[two] = -wheight + radius[two];
								
								//Check for collision with the west wall
								if (xPosition[one] - radius[one] < -wlength + wwidth)
									xPosition[one] = -wlength + wwidth + radius[one];
								if (xPosition[two] - radius[two] < -wlength + wwidth)
									xPosition[two] = -wlength + wwidth + radius[two];
								//Check for collision with the east wall
								if (xPosition[one] + radius[one] > wlength - wwidth)
									xPosition[one] = wlength - wwidth -radius[one];
								if (xPosition[two] + radius[two] > wlength - wwidth)
									xPosition[two] = wlength - wwidth -radius[two];

								//if the distance between the balls are larger than the sum of their radius
								if(sqrtf(pow(xPosition[one]-xPosition[two], 2) + pow(yPosition[one]-yPosition[two], 2)) > radius[one] + radius[two])
								{
									xPosition[one] = tempXPos[0];
									yPosition[one] = tempYPos[0];

									xPosition[two] = tempXPos[1];
									yPosition[two] = tempYPos[1];
									break;
								}

							}
							//the collision can now be calculated since the distance between them is ok

							//initialize some parameters	
							float totMass = mass[one]+mass[two];
							posOne = vec2(xPosition[one], yPosition[one]);
							posTwo = vec2(xPosition[two], yPosition[two]);
							norm = posOne - posTwo;

							//calculate both normal plane and the collision plane
							normalPlane = normalize(norm);
							vec2 collisionPlane = vec2(-normalPlane[1], normalPlane[0]);

							//Calculate prior velocities relative the the collision plane and normal
							vec2 v1 = vec2(xVel[one], yVel[one]);
							vec2 v2 = vec2(xVel[two], yVel[two]);
							float n_vel1 = dot(normalPlane, v1);
							float c_vel1 = dot(collisionPlane, v1);
							float n_vel2 = dot(normalPlane, v2);
							float c_vel2 = dot(collisionPlane, v2);

							// Calculate the scalar velocities of each object after the collision.
							float n_vel1_after = ((n_vel1 * (mass[one] - mass[two])) + (2 * mass[two] * n_vel2)) / (totMass);
							float n_vel2_after = ((n_vel2 * (mass[two] - mass[one])) + (2 * mass[one] * n_vel1)) / (totMass);
							float velObject1Tangent_After = c_vel1;
							float velObject2Tangent_After = c_vel2;

							// Convert the scalars to vectors by multiplying by the normalised plane vectors.
							vec2 vec_n_vel2_after = n_vel2_after * normalPlane;
							vec2 vec_c_vel2 = c_vel2 * collisionPlane;
							vec2 vec_n_vel1_after = n_vel1_after * normalPlane;
							vec2 vec_c_vel1 = c_vel1 * collisionPlane;

							// Combine the vectors back into a single vector in world space.
							vec2 vel1_after = vec_n_vel1_after + vec_c_vel1;
							vec2 vel2_after = vec_n_vel2_after + vec_c_vel2;

							xVel[one] = vel1_after[0]*0.9;
							xVel[two] = vel2_after[0]*0.9;

							yVel[one] = vel1_after[1]*0.9;
							yVel[two] = vel2_after[1]*0.9;
						}
					}
				}
			}while(collision == true);
		}

		//draw the ballz -> since all balls, collided or not, have accurate values. 
		for(int index = 0; index < n_ballz; index++){
			int i = 189*index;
			float ang = 0.0f;
			float step = 0.1f; //same decimals for even size of g_vertex_buffer_data[int]
			float pi2 = 3.1f;

			//cout << index << " x: " << xVel[index]  << " y: " << yPosition[index] << endl;
      		
      		//draw each ball: consist of vertices calculated with sin and cos.
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

		//bind the buffers for the balls: vertices and the colors
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

		//buffers for the walls and ground
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

		// Draw the walls
		glDrawArrays(GL_TRIANGLES, 0, 3*6); // 3 indices starting at 0 -> 1 triangle

		//buffers for the rope
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
		//draw the ropes
		glDrawArrays(GL_LINES, 0, 2*n_ballz);
		//activate attribute(1):colors
		glEnableVertexAttribArray(1); 
    
		//bind the balls
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
		//bind the colors: attribute(1)->add color to the object
		glBindBuffer(GL_ARRAY_BUFFER, colorBallBuffer); 
		glVertexAttribPointer(
								1,                  // attribute 1. No particular reason for 0, but must match the layout in the shader.
								3,                  // size
								GL_FLOAT,           // type
								GL_FALSE,           // normalized?
								0,                  // stride
								(void*)0            // array buffer offset
								);

		// Draw the balls!
		for(int w = 0; w < n_ballz; w++)
			glDrawArrays(GL_TRIANGLE_FAN, w*63, 63); // 3 indices starting at 0 -> 1 triangle
	

		//Disable attributes
		glDisableVertexAttribArray(1); 
		glDisableVertexAttribArray(0);
			
		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();
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
