/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cn.edu.hit.dnmfilterindel.file;

/**
 *
 * @author Yongzhuang Liu
 */
public class Trio {

    private String familyID;
    private Individual father;
    private Individual mother;
    private Individual offspring;

    public Trio(String familyID, Individual father, Individual mother, Individual offspring) {
        this.familyID = familyID;
        this.father = father;
        this.mother = mother;
        this.offspring = offspring;
    }

    public String getFamilyID() {
        return familyID;
    }

    public Individual getFather() {
        return father;
    }

    public Individual getMother() {
        return mother;
    }

    public Individual getOffspring() {
        return offspring;
    }

    public void setMother(Individual mother) {
    	this.mother = mother;
    }
    
    public void setOffspring(Individual offspring) {
    	this.offspring = offspring;
    }
    
    public void setFather(Individual father) {
    	this.father = father;
    }
}
